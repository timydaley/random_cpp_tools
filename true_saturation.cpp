/*    true_saturation:
 *
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith
 *                       Timothy Daley
 *
 *    Authors: Andrew D. Smith and Timothy Daley
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <smithlab_os.hpp>


#include <fstream>
#include <numeric>
#include <vector>
#include <iomanip>
#include <queue>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <tr1/unordered_map>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>


using std::string;
using std::min;
using std::vector;
using std::endl;
using std::cerr;
using std::max;
using std::ifstream;

using std::setw;
using std::fixed;
using std::setprecision;
using std::priority_queue;
using std::tr1::unordered_map;


//////////////////////////////////////////////////////////////////////
// Data imputation
/////////////////////////////////////////////////////////////////////

static bool
update_pe_duplicate_counts_hist(const GenomicRegion &curr_gr,
                                const GenomicRegion &prev_gr,
                                vector<double> &counts_hist,
                                size_t &current_count){
    // check if reads are sorted
    if (curr_gr.same_chrom(prev_gr) &&
        curr_gr.get_start() < prev_gr.get_start()
        && curr_gr.get_end() < prev_gr.get_end()){
        return false;
    }
    if (!curr_gr.same_chrom(prev_gr) ||
        curr_gr.get_start() != prev_gr.get_start() ||
        curr_gr.get_end() != prev_gr.get_end())
        // next read is new, update counts_hist to include current_count
    {
        // histogram is too small, resize
        if(counts_hist.size() < current_count + 1)
            counts_hist.resize(current_count + 1, 0.0);
        ++counts_hist[current_count];
        current_count = 1;
    }
    else // next read is same, update current_count
        ++current_count;
    
    return true;
}


/*
 * This code is used to deal with read data in BAM format.
 */
#ifdef HAVE_SAMTOOLS
// switching dependency on bamtools to samtools
#include <SAM.hpp>


const string mapper = "general";


static void
update_se_duplicate_counts_hist(const MappedRead &curr_mr,
                                const MappedRead &prev_mr,
                                const string input_file_name,
                                vector<double> &counts_hist,
                                size_t &current_count){
    // check if reads are sorted
    if (curr_mr.r.same_chrom(prev_mr.r) &&
        curr_mr.r.get_start() < prev_mr.r.get_start())
        throw SMITHLABException("locations unsorted in: " + input_file_name);
    
    if (!curr_mr.r.same_chrom(prev_mr.r) ||
        curr_mr.r.get_start() != prev_mr.r.get_start())
        // next read is new, update counts_hist to include current_count
    {
        // histogram is too small, resize
        if(counts_hist.size() < current_count + 1)
            counts_hist.resize(current_count + 1, 0.0);
        ++counts_hist[current_count];
        current_count = 1;
    }
    else // next read is same, update current_count
        ++current_count;
}


static size_t
load_counts_BAM_se(const string &input_file_name, vector<double> &counts_hist) {
    
    SAMReader sam_reader(input_file_name, mapper);
    if(!(sam_reader.is_good()))
        throw SMITHLABException("problem opening input file " + input_file_name);
    
    SAMRecord samr;
    sam_reader >> samr;
    size_t n_reads = 1;
    // resize vals_hist, make sure it starts out empty
    counts_hist.clear();
    counts_hist.resize(2, 0.0);
    size_t current_count = 1;
    
    MappedRead prev_mr, curr_mr;
    prev_mr = samr.mr;
    
    while (sam_reader >> samr, sam_reader.is_good())
    {
        if(samr.is_primary && samr.is_mapped)
            // only convert mapped and primary reads
        {
            // ignore unmapped reads & secondary alignments
            if(!(samr.is_mapping_paired) ||
               (samr.is_mapping_paired && samr.is_Trich)){
                //only count unpaired reads or the left mate of paired reads
                
                curr_mr = samr.mr;
                update_se_duplicate_counts_hist(curr_mr, prev_mr, input_file_name,
                                                counts_hist, current_count);
                
                // update number of reads and prev read
                ++n_reads;
                prev_mr = samr.mr;
            }
        }
    }
    
    // to account for the last read compared to the one before it.
    if(counts_hist.size() < current_count + 1)
        counts_hist.resize(current_count + 1, 0.0);
    ++counts_hist[current_count];
    
    return n_reads;
}

/********Below are functions for merging pair-end reads********/

static void
merge_mates(const size_t suffix_len, const size_t range,
            const GenomicRegion &one, const GenomicRegion &two,
            GenomicRegion &merged, int &len) {
    
    assert(one.same_chrom(two));
    const size_t read_start = min(one.get_start(), two.get_start());
    const size_t read_end = max(one.get_end(), two.get_end());
    
    len = read_end - read_start;
    
    if(len < 0){
        cerr << one << endl;
        cerr << two << endl;
        throw SMITHLABException("error merging reads");
    }
    
    merged = one;
    merged.set_start(read_start);
    merged.set_end(read_end);
    merged.set_score(one.get_score() + two.get_score());
    
    const string name(one.get_name());
    merged.set_name("FRAG:" + name.substr(0, name.size() - suffix_len));
    
}

// check if reads have same name & chrom
inline static bool
same_read(const size_t suffix_len, 
	  const MappedRead &a, const MappedRead &b) {
  const string sa(a.r.get_name());
  const string sb(b.r.get_name());
  return (std::equal(sa.begin(), sa.end() - suffix_len, sb.begin())
	  && a.r.same_chrom(b.r));
}


/////comparison function for priority queue/////////////////

/**************** FOR CLARITY BELOW WHEN COMPARING READS *************/
static inline bool
chrom_greater(const GenomicRegion &a, const GenomicRegion &b) {
    return a.get_chrom() > b.get_chrom();
}
static inline bool
same_start(const GenomicRegion &a, const GenomicRegion &b) {
    return a.get_start() == b.get_start();
}
static inline bool
start_greater(const GenomicRegion &a, const GenomicRegion &b) {
    return a.get_start() > b.get_start();
}
static inline bool
end_greater(const GenomicRegion &a, const GenomicRegion &b) {
    return a.get_end() > b.get_end();
}
/******************************************************************************/


struct GenomicRegionOrderChecker {
    bool operator()(const GenomicRegion &prev, const GenomicRegion &gr) const {
        return start_check(prev, gr);
    }
    static bool
    is_ready(const priority_queue<GenomicRegion, vector<GenomicRegion>, GenomicRegionOrderChecker> &pq,
             const GenomicRegion &gr, const size_t max_width) {
        return !pq.top().same_chrom(gr) || pq.top().get_end() + max_width < gr.get_start();
    }
    static bool
    start_check(const GenomicRegion &prev, const GenomicRegion &gr) {
        return (chrom_greater(prev, gr)
                || (prev.same_chrom(gr) && start_greater(prev, gr))
                || (prev.same_chrom(gr) && same_start(prev, gr) && end_greater(prev, gr)));
    }
};



static void empty_pq(GenomicRegion &curr_gr, GenomicRegion &prev_gr,
                     size_t &current_count, vector<double> &counts_hist,
                     priority_queue<GenomicRegion, vector<GenomicRegion>,
                                    GenomicRegionOrderChecker> &read_pq,
                     const string &input_file_name ){
    
    curr_gr = read_pq.top();
    //	       cerr << "outputting from queue : " << read_pq.top() << endl;
    read_pq.pop();
    
    //update counts hist
    bool UPDATE_SUCCESS
    = update_pe_duplicate_counts_hist(curr_gr, prev_gr, counts_hist,
                                      current_count);
    if(!UPDATE_SUCCESS){
        //cerr << "prev = " << prev_gr << endl;
        //cerr << "curr = " << curr_gr << endl;
        //cerr << "priority queue : " << endl;
        //while(	 !(read_pq.empty()) ){
          //  cerr << read_pq.top() << endl;
          //  read_pq.pop();
        //}
        throw SMITHLABException("reads unsorted in " + input_file_name);
    }
    prev_gr = curr_gr;
}



static size_t
load_counts_BAM_pe(const bool VERBOSE,
                   const string &input_file_name,
                   const size_t MAX_SEGMENT_LENGTH,
                   const size_t MAX_READS_TO_HOLD,
                   vector<double> &counts_hist) {
    
    SAMReader sam_reader(input_file_name, mapper);
    // check sam_reader
    if(!(sam_reader.is_good()))
        throw SMITHLABException("problem opening input file " + input_file_name);
    
    SAMRecord samr;
    size_t n_reads = 0;
    // resize vals_hist, make sure it starts out empty
    counts_hist.clear();
    counts_hist.resize(2, 0.0);
    size_t current_count = 1;
    size_t suffix_len = 0;
    size_t n_merged = 0;
    size_t n_unpaired = 0;
    
    GenomicRegion curr_gr, prev_gr;
    
    std::priority_queue<GenomicRegion, vector<GenomicRegion>,
    GenomicRegionOrderChecker> read_pq;
    
    unordered_map<string, SAMRecord> dangling_mates;
    
    while ((sam_reader >> samr, sam_reader.is_good()))
    {
      if(samr.is_primary && samr.is_mapped){
	// only convert mapped and primary reads
	if (samr.is_mapping_paired){
	  const string read_name
	    = samr.mr.r.get_name().substr(0, samr.mr.r.get_name().size() - suffix_len);
                
	  if (dangling_mates.find(read_name) != dangling_mates.end()){
	    // other end is in dangling mates, merge the two mates
	    if(same_read(suffix_len, samr.mr, dangling_mates[read_name].mr)){
	      if (samr.is_Trich) std::swap(samr, dangling_mates[read_name]);
                    
	      GenomicRegion merged;
	      int len = 0;
	      merge_mates(suffix_len, MAX_SEGMENT_LENGTH,
			  dangling_mates[read_name].mr.r, samr.mr.r, merged, len);
	    // merge success!
	      if (len >= 0 && len <= static_cast<int>(MAX_SEGMENT_LENGTH)){
	      // first iteration
		if(n_reads == 0){
		  prev_gr = merged;
		  ++n_reads;
		  ++n_merged;
		}
		else{
		  ++n_reads;
		  ++n_merged;
		  read_pq.push(merged);
                            
		  if(!(read_pq.empty()) &&
		     GenomicRegionOrderChecker::is_ready(read_pq, merged, MAX_SEGMENT_LENGTH)) {
		  //begin emptying priority queue
		    while(!(read_pq.empty()) &&
			  GenomicRegionOrderChecker::is_ready(read_pq, merged,
							      MAX_SEGMENT_LENGTH) ){
		      empty_pq(curr_gr, prev_gr, current_count,
			       counts_hist, read_pq, input_file_name);
		    }//end while loop
		  }//end statement for emptying priority queue
		}
		dangling_mates.erase(read_name);
	     

	      }//end if statement for if merge is successful
	      else{
		cerr << "problem with read " << read_name << endl;

		throw SMITHLABException("merge unsuccessful");
	      }
	    }
	    else{
	      // problem mergin reads from "different" chrom like chr1 & chr1_gl000191_random
	      // flagged as proper pair, but not
	      read_pq.push(samr.mr.r);
	      read_pq.push(dangling_mates[read_name].mr.r);
	    }
	  }//end if statement for if read is in dangling mates
	  else{	// other end not in dangling mates, add this read to dangling mates.
	    dangling_mates[read_name] = samr;
	  }
	}
	else{ // read is unpaired, put in queue
	  if(n_reads == 0){
	    ++n_reads;
	    ++n_unpaired;
	    prev_gr = samr.mr.r;
	  }
	  else{ 
	    ++n_reads;
	    ++n_unpaired;
	    read_pq.push(samr.mr.r);
 
	    if(!(read_pq.empty()) &&
	       GenomicRegionOrderChecker::is_ready(read_pq, samr.mr.r, MAX_SEGMENT_LENGTH)) {
	     
	      while(!(read_pq.empty()) &&
		    GenomicRegionOrderChecker::is_ready(read_pq, samr.mr.r,
							MAX_SEGMENT_LENGTH) ){
                              
		empty_pq(curr_gr, prev_gr, current_count,
			 counts_hist, read_pq, input_file_name);
	      }
                    
	    }
                
	  }

	} 
            
	// dangling mates is too large, flush dangling_mates of reads
	// on different chroms and too far away
	if (dangling_mates.size() > MAX_READS_TO_HOLD){
	  if(VERBOSE)
	    cerr << "dangling mates too large, emptying" << endl;
                
	  unordered_map<string, SAMRecord> tmp;
	  for (unordered_map<string, SAMRecord>::iterator
		 itr = dangling_mates.begin();
	       itr != dangling_mates.end(); ++itr)
	    if (itr->second.mr.r.get_chrom() != samr.mr.r.get_chrom()
		|| (itr->second.mr.r.get_chrom() == samr.mr.r.get_chrom()
		    && itr->second.mr.r.get_end() + MAX_SEGMENT_LENGTH <
		    samr.mr.r.get_start())) 
	      {
		if(itr->second.seg_len >= 0)
		  read_pq.push(itr->second.mr.r);
	      }
	    else
	      tmp[itr->first] = itr->second;
                
	  std::swap(tmp, dangling_mates);
	}
            
      }
    }
    
    // empty dangling mates of any excess reads
    while (!dangling_mates.empty()) {
      read_pq.push(dangling_mates.begin()->second.mr.r);
      dangling_mates.erase(dangling_mates.begin());
    }
  
    //final iteration
    while(!read_pq.empty()){
      empty_pq(curr_gr, prev_gr, current_count,
	       counts_hist, read_pq, input_file_name);
    }
    
    if(counts_hist.size() < current_count + 1)
      counts_hist.resize(current_count + 1, 0.0);

    ++counts_hist[current_count];

    assert((read_pq.empty()));
    
    return n_reads;
}

#endif


/*
 this code is for BED file input
 */

static void
update_se_duplicate_counts_hist(const GenomicRegion &curr_gr,
                                const GenomicRegion &prev_gr,
                                const string input_file_name,
                                vector<double> &counts_hist,
                                size_t &current_count){
    // check if reads are sorted
    if (curr_gr.same_chrom(prev_gr) &&
        curr_gr.get_start() < prev_gr.get_start())
        throw SMITHLABException("locations unsorted in: " + input_file_name);
    
    if (!curr_gr.same_chrom(prev_gr) ||
        curr_gr.get_start() != prev_gr.get_start())
        // next read is new, update counts_hist to include current_count
    {
        // histogram is too small, resize
        if(counts_hist.size() < current_count + 1)
            counts_hist.resize(current_count + 1, 0.0);
        ++counts_hist[current_count];
        current_count = 1;
    }
    else // next read is same, update current_count
        ++current_count;
}

static size_t
load_counts_BED_se(const string input_file_name, vector<double> &counts_hist) {
    // resize vals_hist
    counts_hist.clear();
    counts_hist.resize(2, 0.0);
    
    std::ifstream in(input_file_name.c_str());
    if (!in)
        throw SMITHLABException("problem opening file: " + input_file_name);
    
    GenomicRegion curr_gr, prev_gr;
    if (!(in >> prev_gr))
        throw SMITHLABException("problem opening file: " + input_file_name);
    
    size_t n_reads = 1;
    size_t current_count = 1;
    
    while (in >> curr_gr) {
        update_se_duplicate_counts_hist(curr_gr, prev_gr, input_file_name,
                                        counts_hist, current_count);
        ++n_reads;
        prev_gr.swap(curr_gr);
    }
    
    
    // to account for the last read compared to the one before it.
    if(counts_hist.size() < current_count + 1)
        counts_hist.resize(current_count + 1, 0.0);
    ++counts_hist[current_count];
    
    return n_reads;
}


static size_t
load_counts_BED_pe(const string input_file_name, vector<double> &counts_hist) {
    
    // resize vals_hist
    counts_hist.clear();
    counts_hist.resize(2, 0.0);
    
    std::ifstream in(input_file_name.c_str());
    if (!in)
        throw SMITHLABException("problem opening file: " + input_file_name);
    
    GenomicRegion curr_gr, prev_gr;
    if (!(in >> prev_gr))
        throw SMITHLABException("problem opening file: " + input_file_name);
    
    size_t n_reads = 1;
    size_t current_count = 1;
    
    //read in file and compare each gr with the one before it
    while (in >> curr_gr) {
        bool UPDATE_SUCCESS =
        update_pe_duplicate_counts_hist(curr_gr, prev_gr,
                                        counts_hist, current_count);
        if(!UPDATE_SUCCESS){
           // cerr << "prev = " << prev_gr << endl;
           // cerr << "curr = " << curr_gr << endl;
            throw SMITHLABException("reads unsorted in " + input_file_name);
        }
        
        ++n_reads;
        prev_gr.swap(curr_gr);
    }
    
    if(counts_hist.size() < current_count + 1)
        counts_hist.resize(current_count + 1, 0.0);
    
    // to account for the last read compared to the one before it.
    ++counts_hist[current_count];
    
    
    return n_reads;
    
    
}

/*
 text file input
 */

static size_t
load_counts(const string input_file_name, vector<double> &counts_hist) {
    
    std::ifstream in(input_file_name.c_str());
    if (!in) // if file doesn't open
        throw SMITHLABException("problem opening file: " + input_file_name);
    
    size_t n_reads = 0;
    while(!in.eof()){
        string buffer;
        getline(in, buffer);
        
        std::istringstream iss(buffer);
        if(iss.good()){
            double val;
            iss >> val;
            if(val > 0) {
                const size_t count = static_cast<size_t>(val);
                // histogram is too small, resize
                if(counts_hist.size() < count + 1)
                    counts_hist.resize(count + 1, 0.0);
                ++counts_hist[count];
                n_reads += count;
            }
            else if(val != 0)
                throw SMITHLABException("problem reading file at line " + (n_reads + 1));
        }
        in.peek();
    }
    in.close();
    
    return n_reads;
}

//returns number of reads from file containing counts histogram
static size_t
load_histogram(const string &filename, vector<double> &counts_hist) {
    
    counts_hist.clear();
    
    std::ifstream in(filename.c_str());
    if (!in) //if file doesn't open
        throw SMITHLABException("could not open histogram: " + filename);
    
    size_t n_reads = 0;
    size_t line_count = 0ul, prev_read_count = 0ul;
    string buffer;
    while (getline(in, buffer)) {
        ++line_count;
        size_t read_count = 0ul;
        double frequency = 0.0;
        std::istringstream is(buffer);
        // error reading input
        if (!(is >> read_count >> frequency))
            throw SMITHLABException("bad histogram line format:\n" +
                                    buffer + "\n(line " + toa(line_count) + ")");
        // histogram is out of order
        if (read_count < prev_read_count)
            throw SMITHLABException("bad line order in file " +
                                    filename + "\n(line " +
                                    toa(line_count) + ")");
        counts_hist.resize(read_count + 1, 0.0);
        counts_hist[read_count] = frequency;
        prev_read_count = read_count;
        n_reads += static_cast<size_t>(read_count*frequency);
    }
    
    return n_reads;
}


double 
log_sum_log_vec(const vector<double> &vals, size_t limit){
  const size_t max_idx = 
  max_element(vals.begin(), vals.begin() + limit) - 
  vals.begin();
  const double max_val = vals[max_idx];
  double sum = 1.0;
  for (size_t i = 0; i < limit; ++i) {
    if (i != max_idx) {
      sum += exp(vals[i] - max_val);
#ifdef DEBUG
      assert(finite(sum)); 
      // abort if the sum is infinte //
#endif
    }
  }
  return(max_val + log(sum));
}


static double
expected_coverage(const vector<double> &vals_hist,
                  const double vals_sum,
                  const double extrapolation_size){
  vector<double> log_probs_vec;
  for(size_t i = 1; i < vals_hist.size(); i++)
    log_probs_vec.push_back(log(vals_hist[i])
			    +log(i) - log(vals_sum)
			    -i*extrapolation_size/vals_sum);
  return(exp(log_sum_log_vec(log_probs_vec, log_probs_vec.size())));
}

int
main(int argc, const char **argv){


  try {
    /* FILES */
    string outfile;
    bool VERBOSE = false;

    double max_extrap = 1e10;
    double step_size = 1e6;
    bool VALS_INPUT = false;
    bool PAIRED_END = false;
    bool HIST_INPUT = false;
    
#ifdef HAVE_SAMTOOLS
    bool BAM_FORMAT_INPUT = false;
    size_t MAX_SEGMENT_LENGTH = 1000000;
#endif   
     
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("true_saturation", "",
			   "*.txt");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("verbose", 'v', "print more information", 
		      false , VERBOSE);
    opt_parse.add_opt("max_extrap", 'x', "maximum extrapolation",
		      false, max_extrap);
    opt_parse.add_opt("step_size", 'p', "step size between estimates",
		      false, step_size);
#ifdef HAVE_SAMTOOLS
    opt_parse.add_opt("bam", 'B', "input is in BAM format", 
		      false, BAM_FORMAT_INPUT);
#endif
    opt_parse.add_opt("pe", 'P', "input is paired end read file",
		      false, PAIRED_END);
    opt_parse.add_opt("vals", 'V', 
		      "input is a text file containing only the observed counts",
		      false, VALS_INPUT);
    opt_parse.add_opt("hist", 'H', 
		      "input is a text file containing the observed histogram",
		      false, HIST_INPUT);


    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string input_file_name = leftover_args.front();
    /**********************************************************************/
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    srand(time(0) + getpid());
    gsl_rng_set(rng, rand());
    
    vector<double> counts_hist;
    size_t n_reads = 0;
    
    // LOAD VALUES
    if(HIST_INPUT){
        if(VERBOSE)
            cerr << "INPUT_HIST" << endl;
        n_reads = load_histogram(input_file_name, counts_hist);
    }
    else if(VALS_INPUT){
        if(VERBOSE)
            cerr << "VALS_INPUT" << endl;
        n_reads = load_counts(input_file_name, counts_hist);
    }
#ifdef HAVE_SAMTOOLS
    else if (BAM_FORMAT_INPUT && PAIRED_END){
        if(VERBOSE)
            cerr << "PAIRED_END_BAM_INPUT" << endl;
        const size_t MAX_READS_TO_HOLD = 100000;
        n_reads = load_counts_BAM_pe(VERBOSE, input_file_name, MAX_SEGMENT_LENGTH,
                                     MAX_READS_TO_HOLD, counts_hist);
    }
    else if(BAM_FORMAT_INPUT){
        if(VERBOSE)
            cerr << "BAM_INPUT" << endl;
        n_reads = load_counts_BAM_se(input_file_name, counts_hist);
    }
#endif
    else if(PAIRED_END){
        if(VERBOSE)
            cerr << "PAIRED_END_BED_INPUT" << endl;
        n_reads = load_counts_BED_pe(input_file_name, counts_hist);
    }
    else{ // default is single end bed file
        if(VERBOSE)
            cerr << "BED_INPUT" << endl;
        n_reads = load_counts_BED_se(input_file_name, counts_hist);
    }

    double vals_sum = 0.0;
    for(size_t i = 0; i < counts_hist.size(); i++)
      vals_sum += i*counts_hist[i];

  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  out << "reads" << '\t' << "expected_saturation" << endl;
  out << 0 << '\t' << 1.0 << endl;
    for(double extrap = step_size;
	extrap <= max_extrap; extrap += step_size){
      out << extrap << "\t" << expected_coverage(counts_hist, vals_sum, extrap) << endl;
    }

  }  
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
