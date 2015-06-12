#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"
#include "MappedRead.hpp"
#include "smithlab_os.hpp"
#include "RNG.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cstdlib>

using std::string;
using std::ios_base;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::getline;
using std::ceil;

static double
get_approx_num_lines(const bool VERBOSE, 
		     const string &filename,
		     const size_t n_samples){
  
  const size_t filesize = get_filesize(filename);

  const size_t increment = 
    ceil(static_cast<double>(filesize)/n_samples);

  if(VERBOSE)
    cerr << "increment             = " << increment << endl;
  
  std::ifstream in(filename.c_str(), ios_base::binary);
  if (!in)
    throw SMITHLABException("cannot open input file " + string(filename));
  
  double sampled_bytes = 0.0;
  size_t taken_samples = 0;
  for (size_t i = 0; i < filesize && in.good(); i += increment) {
    Runif ran(time(0) + getpid());
    int ran_bytes = ran.runif(i, std::min(i + increment, filesize));
    in.seekg(ran_bytes, ios_base::beg);
    in.ignore(10000ul, '\n');
    string sample_string;
    getline(in, sample_string);
    if (in.good()){
      sampled_bytes += (sample_string.size() + 1.0)*sizeof(char);
      ++taken_samples;
    }
  }
  return filesize/(sampled_bytes/taken_samples);
}


static size_t
num_bytes_to_next_sample(const gsl_rng *rng,
			 const double prob,
			 const size_t bytes_per_line){
  const unsigned int n_lines = gsl_ran_geometric(rng, prob);
  return static_cast<size_t>(n_lines)*bytes_per_line;
}


int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    size_t sample_size = 1000000;
    size_t n_samples = 1000;
    const size_t buffer_size = 10000;


    bool VERBOSE = false;


    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("sample_reads_noreplace", "sample reads w/out replacement, where each read is a line",
			   "<bed-format-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("sample_size", 's', "size of sample",
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

    const size_t filesize = get_filesize(input_file_name);
    const double approx_num_lines = get_approx_num_lines(VERBOSE, input_file_name, n_samples);
    const double bytes_per_line = filesize/approx_num_lines;
    const double prob = 
      static_cast<double>(sample_size)/ceil(approx_num_lines);

    if(VERBOSE) {
      cerr << "input file size       = " << filesize << endl
	   << "approx bytes per line = " << bytes_per_line << endl
	   << "approx num lines      = " << approx_num_lines << endl
	   << "prob                  = " << prob << endl;
    }

    std::ifstream in(input_file_name.c_str());
    if (!in)
      throw "problem opening file: " + input_file_name;

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    size_t n_out = 0ul;
    size_t bytes_gone = num_bytes_to_next_sample(rng, prob, bytes_per_line);
    while(bytes_gone < filesize && !in.eof() && in.good()){
      // skip to next chosen line
      in.seekg(bytes_gone, ios_base::beg);
      in.ignore(buffer_size, '\n');
      // get line to output
      string out_string;
      getline(in, out_string);
      if(in.good())
	out << out_string << endl;
      ++n_out;
      // update next chosen line
      bytes_gone += num_bytes_to_next_sample(rng, prob, bytes_per_line);
    }
    if(VERBOSE)
      cerr << "num sampled lines     = " << n_out << endl;

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
