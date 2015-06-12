#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"
#include "MappedRead.hpp"

#include <gsl/gsl_rng.h>
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


static size_t
load_values(const string input_file_name, vector<double> &values) {

  std::ifstream in(input_file_name.c_str());
  if (!in)
    throw SMITHLABException("problem opening file: " + input_file_name);

  vector<double> full_values;
  size_t n_reads = 0;
  static const size_t buffer_size = 10000; // Magic!
  while(!in.eof()){
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    double val = atof(buffer);
    if(val > 0.0)
    full_values.push_back(val);
    if(full_values.back() < 0.0){
      cerr << "INVALID INPUT\t" << buffer << endl;
      throw SMITHLABException("ERROR IN INPUT");
    }
    ++n_reads;
    in.peek();
  }
  in.close();
  if(full_values.back() == 0)
    full_values.pop_back();

  values.swap(full_values);
  return n_reads;
}

static void
sample_vals_replace(const gsl_rng *rng,
                    const vector<double> &full_values,
                    const size_t sample_size,
                    vector<double> &sample_vals){

  vector<unsigned int> sample_counts(full_values.size(), 0);
  gsl_ran_multinomial(rng, full_values.size(),
                      static_cast<unsigned int>(sample_size),
                      &full_values.front(), &sample_counts.front());

  sample_vals.clear();
  for(size_t i = 0; i < sample_counts.size(); i++)
    if(sample_counts[i] > 0)
      sample_vals.push_back(static_cast<double>(sample_counts[i]));

}



static void
sample_vals_noreplace(const gsl_rng *rng,
		    const vector<double> &full_values,
		    const size_t sample_size,
		    vector<double> &sample_vals){

  vector<size_t> orig;
  for (size_t i = 0; i < full_values.size(); ++i){
    for(size_t j = 0; j < full_values[i]; j++)
      orig.push_back(i);
  }

  vector<size_t> sample_orig(sample_size, 0);
  gsl_ran_choose(rng, (size_t *)&sample_orig.front(), sample_size, 
		 (size_t *)&orig.front(), orig.size(), sizeof(size_t));

  vector<double> sample_counts;
  double count  = 1;
  for(size_t i = 1; i < sample_orig.size(); i++){
    if(sample_orig[i] != sample_orig[i - 1]){
      sample_counts.push_back(count);
      count = 1;
    }
    else
      count++;
  }
  sample_counts.push_back(count);

  sample_vals.swap(sample_counts);
}

int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    size_t sample_size = 1000000;

    bool VERBOSE = false;
    bool REPLACE = false;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("counts.txt", "",
			   "<count-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("sample_size", 's', "size of sample",
                      false, sample_size);
    opt_parse.add_opt("REPLACE", 'R', "include if sampling is to be done with replacement",
		      false, REPLACE);
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
    vector<double> full_values;
    load_values(input_file_name, full_values);

    if(VERBOSE){
      const size_t max_observed_count = 
	static_cast<size_t>(*std::max_element(full_values.begin(), 
					      full_values.end()));
    
    // BUILD THE HISTOGRAM
      vector<double> counts_hist(max_observed_count + 1, 0.0);
      for (size_t i = 0; i < full_values.size(); ++i)
	++counts_hist[static_cast<size_t>(full_values[i])];

      cerr << "TOTAL READS     = " << accumulate(full_values.begin(),
						 full_values.end(), 0.0) << endl
	   << "DISTINCT READS  = " << full_values.size() << endl
	   << "MAX COUNT       = " << max_observed_count << endl
	   << "COUNTS OF 1     = " << counts_hist[1] << endl;

      cerr << "OBSERVED COUNTS (" << counts_hist.size() << ")" << endl;
      for (size_t i = 0; i < counts_hist.size(); i++)
	if (counts_hist[i] > 0)
	  cerr << i << '\t' << counts_hist[i] << endl;
      cerr << endl;
    
    }

  //setup rng
    srand(time(0) + getpid());
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, rand()); 
    vector<double> sample_vals;
    if(REPLACE)
      sample_vals_replace(rng, full_values, sample_size, sample_vals);
    else
      sample_vals_noreplace(rng, full_values, sample_size, sample_vals);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    for(size_t i = 0; i < sample_vals.size(); i++)
      if(sample_vals[i] > 0)
	out << sample_vals[i] << endl;

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
