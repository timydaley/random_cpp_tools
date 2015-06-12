#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"
#include "MappedRead.hpp"
#include "RNG.hpp"

#include <gsl/gsl_randist.h>

#include <fstream>
#include <iomanip>
#include <numeric>
#include <vector>

using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;
using std::pair;
using std::make_pair;
using std::sort;

using std::max;
using std::setw;
using std::fabs;
using std::ceil;
using std::greater;
using std::numeric_limits;


// simple sampling w/ replacement since we know # reads
static void
sample_indexes_replace(const size_t n_reads,
		       const size_t seed,
		       const size_t sample_size,
		       vector<size_t> &sample_indexes){
  sample_indexes.clear();
  Runif ran(seed);
  const size_t start_index = 0;

  for(size_t i = 0; i < sample_size; i++)
    sample_indexes.push_back(ran.runif(start_index, n_reads));

  sort(sample_indexes.begin(), sample_indexes.end());
}


// resevoir sampling to sample w/out replacement
static void
sample_indexes_noreplace(const size_t n_reads,
			 const size_t seed,
			 const size_t sample_size,
			 vector<size_t> &sample_indexes){
  sample_indexes.clear();
  Runif ran(seed);
  const size_t start_index = 0;

  for(size_t i = 0; i < sample_size; i++)
    sample_indexes.push_back(i);

  for(size_t i = sample_size; i < n_reads; i++){
    size_t ran_index = ran.runif(start_index, i);

    if(ran_index < sample_indexes.size())
      sample_indexes[ran_index] = i;
  }

  sort(sample_indexes.begin(), sample_indexes.end());
}

int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    size_t sample_size = 1000000;
    size_t seed = time(0) + getpid();

    bool VERBOSE = false;
    bool REPLACE = false;
    bool PAIRED_END = false;


    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("sample_fastq", "sample reads w/out replacement, where each read is 4 lines",
			   "<fastq-format-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("ran_seed", 's', "seed for random number seed generation (default: time + getpid)", 
		      false, seed);
    opt_parse.add_opt("replace", 'R', "sample reads with replacement (default w/out replacement)",
		      false, REPLACE);
    opt_parse.add_opt("pe", 'P', "paired end data",
		      false, PAIRED_END);
    opt_parse.add_opt("sample_size", 'n', "size of sample",
                      false, sample_size);
    opt_parse.add_opt("verbose", 'v', "print more run information", 
		      false , VERBOSE);
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
    
    // Setup the random number generator
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, time(0) + getpid());
    
    // count lines of input file
    std::ifstream count_lines(input_file_name.c_str());
    if (!count_lines)
      throw SMITHLABException("problem opening file: " + input_file_name);
    size_t num_lines = 0;
    string trash;
    while(getline(count_lines, trash))
      ++num_lines;
    count_lines.close();

    if((num_lines % 4) != 0)
      throw SMITHLABException("number of lines is not divisible by 4");

    const size_t total_reads = num_lines/4;

    if(VERBOSE)
      cerr << "total reads = " << total_reads << endl;

    vector<size_t> sample_indexes;
    if(REPLACE)
      sample_indexes_replace(total_reads, seed, sample_size, sample_indexes);

    else
      sample_indexes_noreplace(total_reads, seed, sample_size, sample_indexes);

    // initialize in & out streams
    std::ifstream fastq_file(input_file_name.c_str());
    if (!fastq_file)
      throw SMITHLABException("problem opening file: " + input_file_name);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    // output sample line numbers
    size_t out_index = 0;
    size_t current_index = 0;
    while(!fastq_file.eof() && out_index < sample_indexes.size()){
      string out_string;
      vector<string> out_read;
      for(size_t i = 0; i < 4; i++){
	getline(fastq_file, out_string);
	out_read.push_back(out_string);
      }
      // test if current index is equal to next sample index
      // output if true, update sample index
      while(current_index == sample_indexes[out_index]){
	for(size_t i = 0; i < out_read.size(); i++)
	  out << out_read[i] << endl;

	++out_index;
      }
      ++current_index;
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
