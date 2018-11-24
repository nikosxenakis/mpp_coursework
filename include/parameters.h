//preprocessor definitions for testing various configurations
//the definitions are made in the Makefile

#ifdef PROD
	#define MAXITER   20000
	#define PRINTFREQ  100
	#define AVERAGE_PIXEL_FOLDER "./data/average_pixel_prod/"
	#define TIME_RESULTS "./data/time_results_prod.tsv"
#endif

#ifdef AVERAGE_PIXEL_TEST
	#define MAXITER   20000
	#define PRINTFREQ  1
	#define AVERAGE_PIXEL_FOLDER "./data/average_pixel_test/"
	#define TIME_RESULTS "./data/time_results_average_pixel_test.tsv"
#endif

#ifdef TIMING_TEST
	#define MAXITER   1500
	#define PRINTFREQ  0
	#define AVERAGE_PIXEL_FOLDER "./data/average_pixel_test/"
	#define TIME_RESULTS "./data/time_results_timing_test.tsv"
#endif

#ifdef TIMING_WITH_INTERVALS_TEST
	#define MAXITER   1500
	#define PRINTFREQ  1
	#define AVERAGE_PIXEL_FOLDER "./data/average_pixel_intervals_test/"
	#define TIME_RESULTS "./data/time_results_timing_intervals_test.tsv"
#endif
