#ifndef __Log_h__
#define __Log_h__

#define LOGGING Log::Instance()

#include <fstream>
#include <ctime>
#include <chrono>
#include <iostream>
#include <iomanip>

#include <boost/filesystem.hpp>
#include <boost/variant.hpp>
#include <boost/range/algorithm.hpp>

#include "hdf5.h"

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "Options.h"
#include "LogMacros.h"

using std::string;

/*
 * Log Singleton
 *
 * Provides global logging functionality
 *
 * Singletons are sometimes frowned-upon - but Logging functionality is one area
 * for which they are considered acceptable.  Doing it this way means the logging
 * class instance doesn't need to be passed around to all and sundry.
 *
 * JR, May 2019
 */


/*
 * Note: I wrote this class when I refactored the SSE code, without reference to the BSE code.
 * While refactoring the BSE code (later) I realised that the logging I had implemented for SSE
 * wasn't sufficient for the logging requirements of the BSE code.  To provide for the logging
 * needs of the BSE code I added new functionality almost as a wrapper around the original, SSE,
 * logging functionality.  Some of the original SSE logging functionality has almost been
 * rendered redundant by the new BSE code, but I have left it here (almost) in its entirety
 * because it may still be useful.  We'll see.
 *
 * When I wrote the SSE logging functionality I provided debugging functionality with it, as
 * well as a set of macros to make debugging and the issuing of warning messages easier.  I
 * also wrote a set of logging macros to make logging easier.  The debug macros are still
 * useful, and I would encourage their use.  As well as writing the BSE logging code, I also
 * wrote the Errors class, which is my attempt at making error handling easier.  Some of the
 * functionality in the Errors class supersedes the DBG_WARN* macros provided here, but the
 * DBG_WARN* macros are still useful in some circumstances (and in fact are still used in
 * various places in the code).  The LOG* macros are somewhat less useful, but are left here
 * in case the original SSE logging functionality (that which underlies the BSE logging
 * functionality) is used in the future (as mentioned above, it could still be useful in
 * some circumstances).
 *
 * JR, September 2019
 * 
 * 
 * 
 * HDF5 File Support
 * =================
 * 
 * Further note: support for HDF5 was added to the existing log functionality, which was written
 * assuming separate files for each of the log file (system parameters, DCOs, SNe, etc.).  That
 * existing functionality has been maintained, and HDF5 added as an option available to the user.
 * HDF5 support in the code has been shoe-horned into the existing functionality, so although
 * there is only a single HDF5 output file (except for detailed output files), each of the groups
 * within the HDF5 file is treated in the logging code as a separate file - because that's how it 
 * was originally written, and we just use the existing infrastructure and fit the HDF5 support 
 * into that framework.  That makes the HDF5 logging code a bit awkward at times, but it works...
 *
 *
 * A brief description of HDF5 files and chunking, in the COMPAS context:
 * (With the caveat that all I know about HDF5 I learned so that I could write this code - it is
 * very likely that there are better ways to do things, so if anyone knows a better way, please
 * either change the code or tell me how to improve it and I'll change it.  Most of what follows
 * is just a brain dump of my reading/research over the past week or so, and it could very well
 * be based on my misunderstanding of what I have read - so if anyone notices something I've
 * misunderstood please let me know so I can improve the code.)
 * 
 *
 * Data in HDF5 files are arranged in groups and datasets.
 *
 *     A COMPAS output file (e.g. BSE_System_Parameters, BSE_RLOF, etc.) maps to an HDF5 group,
 *     where the group name is the name of the COMPAS output file.
 *
 *     A column in a COMPAS output file (e.g. SEED, Mass(1), Radius(2), etc.) maps to an HDF5 dataset,
 *     where the dataset name is the column heading string.
 *
 *     COMPAS column datatype strings are encoded in the dataset meta-details (dataset.dtype).
 *     COMPAS column units strings are attached to HDF5 datasets as attributes.
 *
 * Each dataset in an HDF5 files is broken into "chunks", where a chunk is defined as a number of dataset
 * entries ("chunks" and "chunking" are HDF5 terms).  In COMPAS, all datasets are 1-d arrays (columns), so
 * a chunk is defined as a number of values in the 1-d array (or column).  Chunking can be enabled or not, 
 * but if chunking is not* enabled a dataset cannot be resized - so if chunking is not enabled the size of 
 * the dataset must be known at the time of creation, and the entire dataset created in one go.  That doesn't 
 * work for COMPAS - even though we know the number of systems being evolved, we don't know the number of 
 * entries we'll have in each of the output log files (and therefore the HDF5 datasests if we're logging to 
 * HDF5 files).  So, we enable chunking.
 * 
 * Chunking can improve, or degrade, performance depending upon how it is implemented - mostly related to
 * the chunk size chosen.
 * 
 * Datasets are stored inside an HDF5 file as a number of chunks - the chunks are not guaranteed (not even
 * likely) to be contiguous in the file or on the storage media (HDD, SSD etc.).  Chunks are mapped/indexed
 * in the HDF5 file using a B-tree, and the size of the B-tree, and therefore the traversal time, depends 
 * directly upon the number of chunks allocated for a dataset - so the access time for a chunk increases as
 * the number of chunks in the dataset increases.  So many small chunks will degrade performance.
 * 
 * Chunks are the unit of IO for HDF5 files - all IO to HDF5 is performed on the basis of chunks.  This means
 * that whenever dataset values are accessed (read or written (i.e. changed)), if the value is not already in
 * memory, the entire chunk containing the value must be read from, or written to, the storage media - even
 * if the dataset value being accessed is the only value in the chunk.  So few large chunks could cause 
 * empty, "wasted", space in the HDF5 files (at the end of datasets) - but they could also adversely affect
 * performance by causing unnecessary IO traffic (although probably not much in the way we access data in COMPAS
 * files).
 * 
 * HDF5 files implement a chunk cache on a per-dataset basis.  The default size of the chunk cache is 1MB, and
 * its maximum size is 32MB.  The purpose of the chunk cache is to reduce storage media IO - even with SSDs,
 * memory access is way faster than storage media access, so the more of the file data that can be kept in
 * memory and maipulated there, the better.  Assuming the datatype of a particular dataset is DOUBLE, and
 * therefore consumes 8 bytes of storage space, at its maximum size the chunk cache for that dataset could hold
 * 4,000,000 values - so a single chuck with 4,000,000 values, two chunks with 2,000,000 values, four with 
 * 1,000,000, and so on.  Caching a single chunk defeats the purpose of the cache, so chunk sizes somewhat less
 * that 4,000,000 would be most appropriate if the chunk cache is to be utilised.  Chunks too big to fit in the
 * cache simply bypass the cache and are read from, or written to, the storage media directly.
 * 
 * However, the chunk cache is really only useful for random access of the dataset.  Most, if not all, of the
 * access in the COMPAS context (including post-creation analyses) is serial - the COMPAS code writes the
 * datasets from top to bottom, and later analyses (generally) read the datasets the same way.  Caching the
 * chunks for serial access just introduces overhead that costs memory (not much, to be sure: up to 32MB per 
 * open dataset), and degrades performance (albeit it a tiny bit).  For that reason I disable the chunk cache
 * in COMPAS - so all IO to/from an HDF5 file in COMPAS is directly to/from the storage media.  (To be clear,
 * post-creation analysis software can disable the cache or not when accessing HDF5 files created by COMPAS - 
 * disabling the cache here does not affect how other software accesses the files post-creation).
 * 
 * So many small chunks is not so good, and neither is just a few very large chunks.  So what's the optimum
 * chunk size?  It depends upon several things, and probably the most important of those are the final size
 * of the dataset and the access pattern.
 * 
 * As mentioned above, we tend to access datasets serially, and generally from to to bottom, so larger chunks 
 * would seem appropriate, but not so large that we generate HDF5 files with lots of unused space.  However, 
 * disk space, even SSD space, is cheap, so trading space against performance is probably a good trade.
 * 
 * Also as mentioned above, we don't know the final size of (most of) the datasets when creating the HDF5 in
 * COMPAS - though we do know the number of systems being generated, which allows us to determine an upper
 * bound for at least some of the datasets (though not for groups such as BSE_RLOF).
 * 
 * One thing we need to keep in mind is that when we create the HDF5 file we write each dataset of a group
 * in the same iteration - this is analogous to writing a single record in (e.g.) a CSV log file (the HDF5
 * group corresponds to the CSV file, and the HDF5 datasets in the group correspond to the columns in the
 * CSV file).  So for each iteration - typically each system evolved (though each timestep for detailed
 * output files) we do as many IOs to the HDF5 file as there are datasets in the group (columns in the file).
 * We are not bound to reading or writing a single chunk at a time - but we are bound to reading or writing
 * an integral multiple of whole chunks at a time.
 * 
 * We want to reduce the number of storage media accesses when writing (or later reading) the HDF5 files, so
 * larger chunk sizes are appropriate, but not so large that we create excessively large HDF5 files that have
 * lots of unused space (bearing in mind the trade-off mentioned above), especially when were evolving just
 * a few systems (rather than millions).
 * 
 * To really optimise IO performance for HDF5 files we'd choose chunk sizes that are close to multiples of
 * storage media block sizes, but I chose not to go down that rathole...
 * 
 * Based on everything written above, and some tests, I've chosen a default chunk size of 100,000 (dataset
 * entries) for all datasets (HDF5_DEFAULT_CHUNK_SIZE in constants.h).  This clearly trades performance against
 * storage space.  For the (current) default logfile record specifications, per-binary logfile space is about
 * 1K bytes, so in the very worst case we will waste some space at the end of an HDF5 output file, but the 
 * performance gain, especially for post-creation analyses, is significant.  Ensuring the number of systems 
 * evolved is an integral multiple of this fixed chunk size will minimise storage space waste.
 * 
 * I have added the program option --hdf5-chunk-size to allow users to specify the chunk size - the option
 * value defaults to HDF5_DEFAULT_CHUNK_SIZE.
 * 
 * I have chosen a minimum chunk size of 1000 (HDF5_MINIMUM_CHUNK_SIZE in constants.h).  If the number of
 * systems being evolved is >= HDF5_MINIMUM_CHUNK_SIZE the chunk size used will be the value of the hdf5-chunk-size
 * program option (either HDF5_DEFAULT_CHUNK_SIZE or a value specified by the user), but if the number of 
 * systems being evolved is < HDF5_MINIMUM_CHUNK_SIZE the chunk size used will be HDF5_MINIMUM_CHUNK_SIZE.
 * This is just so we don't waste too much storage space when running small tests - and if they are that small
 * performance is probably not going to be much of an issue, so no real trade-off against storage space.  
 * Detailed output files will use a chunk size of HDF5_MINIMUM_CHUNK_SIZE on the basis that detailed output
 * files generally don't have many thousands of records.
 * 
 * IO to HDF5 files is buffered in COMPAS - we buffer a number of chunks for each open dataset and write the
 * buffer to the file when the buffer fills (or a partial buffer upon file close if the buffer is not full).
 * This IO buffering is not HDF5 or filesystem buffering - this is a COMPAS-internal implementation to improve
 * performance.  We could do the same for the other logfile types one day, but I suspect HDF5 is mostly
 * what people will use - and since each record written to a logfile of type other than HDF5 includes values
 * for all columns in the file, there are fewer IO operations to logfiles of type other than HDF5 (1 per system
 * being evolved) and so IO to logfiles of type other than HDF5 has a far less significant impact on performance
 * than does IO to HDF5 logfiles (where column (dataset) data is written individually).
 * 
 * The default HDF5 IO buffer size is defined in constants.h - HDF5_DEFAULT_IO_BUFFER_SIZE - and I have set it
 * to just 1 - so by default no buffering happens.  I have added the program option --hdf5-buffer-size to allow 
 * users to specify the buffer size - the option value defaults to HDF5_DEFAULT_IO_BUFFER_SIZE.  The HDF5 IO 
 * buffer size is specified as a number of chunks.  Users should increase the buffer size for better performance
 * if memory space allows it.
 * 
 * Users should bear in mind that the combination of HDF5 chunk size and HDF5 IO buffer size affect performance,
 * storage space, and memory usage - so they may need to experiment to find a balance that suits their needs.
 * 
 * 
 * String values stored in HDF5 files
 * ==================================
 *
 * COMPAS writes string data to HDF5 output files as C-type strings.  Python interprets C-type
 * strings in HDF5 files as byte arrays - regardless of the specified datatype when written (COMPAS
 * writes the strings as ASCII data (as can be seen with h5dump), but Python ignores that).  Note that
 * this affects the values in datasets (and attributes) only, not the dataset names (or group names,
 * attribute names, etc.).
 * 
 * The only real impact of this is that if the byte array is printed by Python, it will be displayed
 * as (e.g.) "b'abcde'" rather than just "abcde".  All operations on the data work as expected - it
 * is just the output that is impacted.  If that's an issue, use .decode('utf-8') to decode the byte
 * array as a Python string variable.  E.g.
 * 
 *     str = h5File[Group][Dataset][0]
 *     str is a byte array and print(str) will display (e.g.) b'abcde'
 * 
 *     but
 * 
 *     str = h5File[Group][Dataset][0].decode('utf-8')
 *     str is a Python string and print(str) will display (e.g.) abcde
 *
 * 
 * JR, January 2021
 * 
 * 
 * 
 * Annotations
 * ===========
 * 
 * We have added functionality to allow users to annotate log files.  The original motivation for 
 * annotation functionality was to enable users to describe the contents of custom grid files, but
 * annotations can be used for any reason.  With the ability to annotate the log files, the users
 * can indicate the origin of various input data - e.g. the user could indicate what IMF was used 
 * to draw initial mass values, or what distribution was used to draw mass ratio (q), etc.  
 * 
 * Annotations are written to log files as columns of data (aka datasets in HDF5 files).
 * Annotations are specified via program options, and so can be specified on the command line as
 * well as in grid files.  Because annotations are written to log files as columns of data, we
 * have provided functionality to provide headers for the columns, as well as the actual data
 * for the columns (the annotations).  Two new program options are provided:
 * 
 * --notes-hdrs 
 * Allows users to specify header strings for annotation columns.  Can only be specified on the
 * command line.  Allows users to specify one or more annotation header strings: this is a vector
 * program option.  Usage is:
 * 
 *      --notes-hdrs headerstr1 headerstr2 headerstr3 ... headerstrN
 * 
 * Note that header strings are separated by a space character.  There is no limit to the number
 * of header strings specified, and the number specified defines the number of annotations allowed.
 *
 * 
 * --notes
 * Allows users to specify annotations.  Can be specified on the command line and in a grid file.
 * Allows users to specify one or more annotation strings: this is a vector program option.  Usage
 * is:
 * 
 *      --notes annotation1 annotation2 annotation3 ... annotationN
 * 
 * Note that annotation strings are separated by a space character.  The number of annotation
 * strings is limited to the number of annotation header strings specified (via the --notes-hdrs
 * program option).  If more annotation strings are specified than header strings, the excess
 * annotation strings will be ignored (and a warning displayed).  Note that when using this notation
 * all annotation strings must be provided: there is no mechanism to allow a default annotation
 * using the fallback method for program options (to the command-line value, then to the COMPAS
 * default) - leaving an annotation string blank would be ambiguous (as to which annotation string
 * had been left blank), and specifying "" as an annotation string would be ambiguous (as to whether
 * the use wanted the annotation string to default, or just be a blank string).
 * 
 * Because this notation could become awkward, and to allow for default annotations, a shorthand
 * notation for vector program options has been provided (see notes in Options.h for details).
 * 
 * Usage using the shorthand notation is:
 * 
 *     --notes-hdrs [headerstr1,headerstr2,headerstr3,...,headerStrN]
 * 
 *     --notes [annotation1,annotation2,annotation3,...,annotationN]
 * 
 * Because the parameters are bounded by the brackets, and delimited by commas (and so are now
 * positional), users can omit specific annotations:
 * 
 *     --notes [,,annotation3,,annotation5]
 * 
 * In the example above, annotations 1, 2, 4, and those beyond annotation 5 have been omitted.
 * Annotations 1, 2 & 4 will default - if they are specified in this manner on a grid line they will
 * default to the correspodning annotation specified on the command line; if they are specified in 
 * this manner on the command line they will default to the COMPAS default annotation (the empty 
 * string).  If the number of annotations expected, as defined by the number of annotation headers 
 * specified via the --notes-hdrs program option is more than 5, then annotations beyond annotation 5
 * (the last annotation actually specified by the user) will default in the same manner as described
 * above.
 * 
 * Note that any spaces in annotation header strings and annotation strings need to be enclosed in
 * quotes, or the shell parser will parse them as separate arguments.  If the logfile type is
 * specified as TXT, then any spaces in annotation header strings and annotation strings need to
 * be enclosed in quotes to avoid the shell parser parsing them as separate arguments, but also
 * need to have enclosing quotes propagated to the logfile, or the spaces will be interpreted as
 * delimiters in the logfile - in this cae, the user will need to add enclosing escaped quote
 * characters ('\"') before adding the enclosing quotes.  e.g.:
 * 
 *     --notes-hdrs [headerstr1,"\"header str 2\"",headerstr3,...,headerStrN]
 * 
 * (Note that this is true of all string program option values - spaces that are to be propagated to
 * TXT logfiles need to be enclosed in escaped quotes).
 * 
 * 
 * JR, October 2021
 * 
 * 
 * 
 * Record Types
 * ============
 * 
 * All standard logfiles, except the switch log files, now have a record type property (column).  The record
 * type property is of type LOGRECORDTYPE, which is a typedef for unsigned int (unsigned int allows up to 
 * 4294967296 different integer record types (per standard log file - that should be plenty...).
 * 
 * The record type property can be used to identify and filter records within a standard log file.  The
 * functionality was introduced primarily to support different types of records in the detailed output files
 * (BSE and SSE), but could be useful for other log files.
 *
 * The current use case for the detailed output log files is to differentiate between records written to the
 * file when the binary and constituent stars are known to be in self-consistent states (that is, the attributes
 * of the binary and constituent stars have been correctly and completely updated), and records written to the
 * file, perhaps mid-timestep, when the binary and/or constituent stars may not be self-consistent.  Since the
 * record type property can take any value in the range 0..4294967295 there is scope to identify many different
 * events or situations, in any of the standard log files.  We may want, for example, to indicate that a detailed 
 * output record was written immediately prior to, or immediately following, a particular event or calculation.  
 * Or we may want to indicate that a supernovae record was written to the supernovae file prior to the SN event, 
 * and another record following the SN event.  Or for the RLOF file we may want to differentiate between records 
 * written pre-MT and post_MT - the possibilities are (almost) limitless.
 * 
 * Each standard log file has its own set of record types, as defined in constants.h (e.g. for the BSE Detailed 
 * Output file, see 'enum class BSE_DETAILED_RECORD_TYPE' in constants.h).  The idea is to use positive integers
 * to specify the record type, then record types can be ORed together to form a record type bitmap which can be
 * checked to determine which record type to print (set by program options).
 *  
 * 
 * JR, August 2022
 * 
 */



/*
 * Format a boost::variant value using a format specification passed as a parameter
 *
 * This is defined as a class for use with boost::apply_visitor().
 * It is only ever used by the Log class, hence the reason it is defined here.
 *
 * This function is applied to a boost::variant value via the boost::apply_visitor() function.
 * The function extracts the underlying primitive value stored in the boost::variant value and formats
 * it using the format string constructed here specifically for the underlying primitive type.
 *
 */
class FormatVariantValue: public boost::static_visitor<string> {
public:
    string operator()(const bool                   v, const string fmtStr) const {
                                                      string fmt = OPTIONS->PrintBoolAsString() ? "%5s" : "%1s";
                                                      string vS  = OPTIONS->PrintBoolAsString() ? (v ? "TRUE " : "FALSE") : (v ? "1" : "0");
                                                      return utils::vFormat(fmt.c_str(), vS.c_str());
                                                   }
    string operator()(const int                    v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const short int              v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const long int               v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const long long int          v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const unsigned int           v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "u"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const unsigned short int     v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "u"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const unsigned long int      v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "u"; return utils::vFormat(fmt.c_str(), v); } // also handles OBJECT_ID (typedef)
    string operator()(const unsigned long long int v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "u"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const float                  v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "e"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const double                 v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "e"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const long double            v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "e"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const string                 v, const string fmtStr) const { string fmt = fmtStr; fmt = "%-" + fmt + "s"; return utils::vFormat(fmt.c_str(), v.c_str()); }
    string operator()(const ERROR                  v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const STELLAR_TYPE           v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const MT_CASE                v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const MT_TRACKING            v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const SN_EVENT               v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const SN_STATE               v, const string fmtStr) const { string fmt = fmtStr; fmt = "%"  + fmt + "d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const std::vector<string>    v, const string fmtStr) const { string fmt = fmtStr; fmt = "%-" + fmt + "s"; return utils::vFormat(fmt.c_str(), v[0].c_str()); }
    string operator()(const std::vector<string>    v, const string fmtStr, const size_t idx) const { string fmt = fmtStr; fmt = "%-" + fmt + "s"; return utils::vFormat(fmt.c_str(), v[idx].c_str()); }
};


/*
 * Format a boost::variant value using a type-based default format specification
 *
 * This is defined as a class for use with boost::apply_visitor().
 * It is only ever used by the Log class, hence the reason it is defined here.
 *
 * This function is applied to a boost::variant value via the boost::apply_visitor() function.
 * The function extracts the underlying primitive value stored in the boost::variant value and formats
 * it using a type-based default format specification.
 *
 */
class FormatVariantValueDefault: public boost::static_visitor<string> {
public:
    string operator()(const bool                   v) const {
                                                      string fmt = OPTIONS->PrintBoolAsString() ? "%5s" : "%1s";
                                                      string vS  = OPTIONS->PrintBoolAsString() ? (v ? "TRUE " : "FALSE") : (v ? "1" : "0");
                                                      return utils::vFormat(fmt.c_str(), vS.c_str());
                                                   }
    string operator()(const int                    v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const short int              v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const long int               v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const long long int          v) const { string fmt = "%28.1d"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const unsigned int           v) const { string fmt = "%14.1u"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const unsigned short int     v) const { string fmt = "%14.1u"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const unsigned long int      v) const { string fmt = "%14.1u"; return utils::vFormat(fmt.c_str(), v); } // also handles OBJECT_ID (typedef)
    string operator()(const unsigned long long int v) const { string fmt = "%28.1u"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const float                  v) const { string fmt = "%16.8e"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const double                 v) const { string fmt = "%16.8e"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const long double            v) const { string fmt = "%16.8e"; return utils::vFormat(fmt.c_str(), v); }
    string operator()(const string                 v) const { string fmt = "%-30s";  return utils::vFormat(fmt.c_str(), v.c_str()); }
    string operator()(const ERROR                  v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const STELLAR_TYPE           v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const MT_CASE                v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const MT_TRACKING            v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const SN_EVENT               v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const SN_STATE               v) const { string fmt = "%14.1d"; return utils::vFormat(fmt.c_str(), static_cast<int>(v)); }
    string operator()(const std::vector<string>    v) const { string fmt = "%-30s"; return utils::vFormat(fmt.c_str(), v[0].c_str()); }
    string operator()(const std::vector<string>    v, const size_t idx) const { string fmt ="%-30s"; return utils::vFormat(fmt.c_str(), v[idx].c_str()); }
};


class Log {

private:

    Log() {                                                                         // constructor - initialise variables
        m_Enabled = false;                                                          // logging disabled initially
        m_HDF5ContainerId = -1;                                                     // no HDF5 container file open initially
        m_Run_Details_H5_File.fileId = -1;                                          // no HDF5 file id for run details file initially
        m_Run_Details_H5_File.groupId = -1;                                         // no HDF5 group id for run details file initially
        m_HDF5DetailedId = -1;                                                      // no HDF5 detailed file open initially
        m_LogBasePath = ".";                                                        // default log file base path
        m_LogContainerName = DEFAULT_OUTPUT_CONTAINER_NAME;                         // default log file container name                        
        m_LogNamePrefix = "";                                                       // default log file name prefix
        m_LogfileType = DEFAULT_LOGFILE_TYPE;                                       // default log file type
        m_LogLevel = 0;                                                             // default log level - log everything
        m_LogClasses = {};                                                          // no default log classes
        m_DbgLevel = 0;                                                             // default debug level - debug everything
        m_DbgClasses = {};                                                          // no default debug classes
        m_DbgToLogfile = false;                                                     // default is not to log debug records to the log file
        m_DbgLogfileId = -1;                                                        // default is not valid
        m_Logfiles.empty();                                                         // default is no log files
        m_OpenStandardLogFileIds = {};                                              // no open COMPAS standard log files

        m_ObjectIdSwitching = -1L;                                                  // object id of Star object switching stellar type - default none
        m_TypeSwitchingFrom = STELLAR_TYPE::NONE;                                   // stellar type from which Star object is switching - default NONE
        m_TypeSwitchingTo   = STELLAR_TYPE::NONE;                                   // stellar type to which Star object is switching - default NONE
        m_PrimarySwitching  = false;                                                // Star switching is primary star of binary - default false

        m_SSESupernovae_DelayedWrite.logRecordType       = 0;                       // delayed log record type for SSE_Supernovae file - initially 0 (set later)
        m_SSESupernovae_DelayedWrite.logRecordString     = "";                      // delayed log record (string) for SSE_Supernovae file - initially empty
        m_SSESupernovae_DelayedWrite.logRecordValues     = {};                      // delayed log record (property values) for SSE_Supernovae file - initially empty
        m_SSESupernovae_DelayedWrite.logRecordProperties = {};                      // SSE Supernovae logfile record properties - initially empty
        m_SSESupernovae_DelayedWrite.logRecordFmtVector  = {};                      // SSE Supernovae logfile format vector - initially empty
        m_SSESupernovae_DelayedWrite.logFileAnnotations  = {};                      // SSE Supernovae annotations vector - initially empty

        m_OptionDetails = {};                                                       // option details retrieved from commandline - initially empty
    };
    Log(Log const&) = delete;                                                       // copy constructor does nothing, and not exposed publicly
    Log& operator = (Log const&) = delete;                                          // operator = does nothing, and not exposed publicly

    // instance variable
    static Log                 *m_Instance;                                         // pointer to the instance


    // member variables
    bool                        m_Enabled;                                          // is logging enabled?

    string                      m_HDF5ContainerName;                                // HDF5 container name
    hid_t                       m_HDF5ContainerId;                                  // HDF5 container id
    hid_t                       m_HDF5DetailedId;                                   // HDF5 detailed output id

    string                      m_LogBasePath;                                      // base path for log files
    string                      m_LogContainerName;                                 // container (directory) name for log files
    string                      m_LogNamePrefix;                                    // prefix for log files

    LOGFILETYPE                 m_LogfileType;                                      // logfile type
    int                         m_LogLevel;                                         // log level
    std::vector <string>        m_LogClasses;                                       // log classes

    int                         m_DbgLevel;                                         // debug level
    std::vector <string>        m_DbgClasses;                                       // debug classes

    bool                        m_DbgToLogfile;                                     // log debug records to log file?
    int                         m_DbgLogfileId;                                     // log file id of file to which debug statements should be written

    bool                        m_ErrToLogfile;                                     // log error records to log file?
    int                         m_ErrLogfileId;                                     // log file id of file to which error statements should be written

    std::vector<OptionDetailsT> m_OptionDetails;                                    // option details retrieved from commandline


    struct h5AttrT {                                                                // attributes of HDF5 files
        hid_t   fileId;                                                             //    - file id
        hid_t   groupId;                                                            //    - group id

        size_t  chunkSize;                                                          //    - chunk size
        size_t  IOBufSize;                                                          //    - IO buffer size

        struct h5DataSetsT {                                                        // attributes of HDF5 datasets
            hid_t                             dataSetId;                            //    - HDF5 dataset id
            hid_t                             h5DataType;                           //    - HDF5 datatype
            TYPENAME                          dataType;                             //    - COMPAS data type
            STRING_QUALIFIER                  stringType;                           //    - Qualifier for string datatype - fixed or variable length
            std::vector<COMPAS_VARIABLE_TYPE> buf;                                  //    - write buffer - for chunking
        };

        std::vector<h5DataSetsT> dataSets;                                          // details of datasets
    };

    struct logfileAttrT {
        bool        active;                                                         // currently logging?  (ie log file open)
        LOGFILE     logfiletype;                                                    // logfile type
        LOGFILETYPE filetype;                                                       // file type
        string      name;                                                           // name of log file
        bool        timestamp;                                                      // time stamp enabled?
        bool        label;                                                          // record labels enabled?

        std::ofstream file;                                                         // file pointer for CSV, TSV, TXT files

        h5AttrT h5File;                                                             // file details for HDF5 files
    };

    std::vector<logfileAttrT> m_Logfiles;                                           // logfiles - in use and not

    h5AttrT m_Run_Details_H5_File;                                                  // HDF5 attributes for run details in HDF5 container

    COMPASUnorderedMap<LOGFILE, LogfileDetailsT> m_OpenStandardLogFileIds;          // currently open standard logfiles: id, filename, property details, field format strings

    // logfile record specifications
    // BSE
    ANY_PROPERTY_VECTOR m_BSE_BE_Binaries_Rec = BSE_BE_BINARIES_REC;                // default specification
    ANY_PROPERTY_VECTOR m_BSE_CEE_Rec         = BSE_COMMON_ENVELOPES_REC;           // default specification
    ANY_PROPERTY_VECTOR m_BSE_DCO_Rec         = BSE_DOUBLE_COMPACT_OBJECTS_REC;     // default specification
    ANY_PROPERTY_VECTOR m_BSE_Detailed_Rec    = BSE_DETAILED_OUTPUT_REC;            // default specification
    ANY_PROPERTY_VECTOR m_BSE_Pulsars_Rec     = BSE_PULSAR_EVOLUTION_REC;           // default specification
    ANY_PROPERTY_VECTOR m_BSE_RLOF_Rec        = BSE_RLOF_PARAMETERS_REC;            // default specification
    ANY_PROPERTY_VECTOR m_BSE_SNE_Rec         = BSE_SUPERNOVAE_REC;                 // default specification
    ANY_PROPERTY_VECTOR m_BSE_Switch_Rec      = BSE_SWITCH_LOG_REC;                 // default specification
    ANY_PROPERTY_VECTOR m_BSE_SysParms_Rec    = BSE_SYSTEM_PARAMETERS_REC;          // default specification

    // SSE
    ANY_PROPERTY_VECTOR m_SSE_Detailed_Rec    = SSE_DETAILED_OUTPUT_REC;            // default specification
    ANY_PROPERTY_VECTOR m_SSE_SNE_Rec         = SSE_SUPERNOVAE_REC;                 // default specification
    ANY_PROPERTY_VECTOR m_SSE_Switch_Rec      = SSE_SWITCH_LOG_REC;                 // default specification
    ANY_PROPERTY_VECTOR m_SSE_SysParms_Rec    = SSE_SYSTEM_PARAMETERS_REC;          // default specification

    // logfile annotation specifications
    //
    // these are just vectors of booleans, each with size equalling the number of annotations defined by the 
    // user-specified option '--notes-hdrs' (OPTIONS->NotesHdrs())
    //
    // the boolean value indicates whether the respective annotation should be recorded in the
    // associated file
    //
    // the default value for each boolean in the vector is FALSE, and the value is only used if the
    // logfile record specification include PROGRAM_OPTION::NOTES.  The defaults are updated in
    // Log::Start(): if PROGRAM_OPTION::NOTES is included in the default record specification for
    // a logfile, the *_Notes defaults are set to TRUE.  This is so Log::UpdateAllLogfileRecordSpecs()
    // has the right defaults when processing any log definitions file.

    // BSE
    std::vector<bool> m_BSE_BE_Binaries_Notes = std::vector<bool>(OPTIONS->NotesHdrs().size(), false);
    std::vector<bool> m_BSE_CEE_Notes         = std::vector<bool>(OPTIONS->NotesHdrs().size(), false);
    std::vector<bool> m_BSE_DCO_Notes         = std::vector<bool>(OPTIONS->NotesHdrs().size(), false);
    std::vector<bool> m_BSE_Detailed_Notes    = std::vector<bool>(OPTIONS->NotesHdrs().size(), false);
    std::vector<bool> m_BSE_Pulsars_Notes     = std::vector<bool>(OPTIONS->NotesHdrs().size(), false);
    std::vector<bool> m_BSE_RLOF_Notes        = std::vector<bool>(OPTIONS->NotesHdrs().size(), false);
    std::vector<bool> m_BSE_SNE_Notes         = std::vector<bool>(OPTIONS->NotesHdrs().size(), false);
    std::vector<bool> m_BSE_Switch_Notes      = std::vector<bool>(OPTIONS->NotesHdrs().size(), false);
    std::vector<bool> m_BSE_SysParms_Notes    = std::vector<bool>(OPTIONS->NotesHdrs().size(), false);

    // SSE
    std::vector<bool> m_SSE_Detailed_Notes    = std::vector<bool>(OPTIONS->NotesHdrs().size(), false);
    std::vector<bool> m_SSE_SNE_Notes         = std::vector<bool>(OPTIONS->NotesHdrs().size(), false);
    std::vector<bool> m_SSE_Switch_Notes      = std::vector<bool>(OPTIONS->NotesHdrs().size(), false);
    std::vector<bool> m_SSE_SysParms_Notes    = std::vector<bool>(OPTIONS->NotesHdrs().size(), false);


    // the following block of variables support the BSE Switch Log file
    
    OBJECT_ID    m_ObjectIdSwitching;                                               // the object id of the Star object switching stellar type
    STELLAR_TYPE m_TypeSwitchingFrom;                                               // the stellar type from which the Star object is switching
    STELLAR_TYPE m_TypeSwitchingTo;                                                 // the stellar type to which the Star object is switching
    bool         m_PrimarySwitching;                                                // flag to indicate whether the primary star of the binary is switching


    // the following struct supports delayed writes to logfiles
    //
    // For now, only delayed writes to the SSE Supernova file is implemented, and only a single record can be delayed.
    // This functionality was introduced specifically to allow queueing a delayed write the to SSE Supernova file - see
    // the discussion in the description of Log::GetStandardLogFileRecordDetails() in Log.cpp.
    // This functionality probably shouldn't be extended to allow queueing/delaying multiple records for later writing
    // (I don't think we need it, it would probably soak up too much memory if over-used, and it might just cause confusion)

    struct delayedWriteDetailsT {                                                   // attributes of delayed writes
        LOGRECORDTYPE                     logRecordType;                            // log record type
        string                            logRecordString;                          // log record to be written to log file in delayed write
        std::vector<COMPAS_VARIABLE_TYPE> logRecordValues;                          // log record property values be written to log file in delayed write
        ANY_PROPERTY_VECTOR               logRecordProperties;                      // logfile record properties
        std::vector<string>               logRecordFmtVector;                       // logfile format vector
        std::vector<bool>                 logFileAnnotations;                       // logfile annotations vector
    };

    delayedWriteDetailsT m_SSESupernovae_DelayedWrite;                              // SSE_Supernovae delayed write details    
    
  
    // the following block of variables support the run details file

    std::ofstream                                      m_RunDetailsFile;            // run details file
    std::chrono::time_point<std::chrono::system_clock> m_WallStart;                 // for run details file
    std::chrono::time_point<std::chrono::system_clock> m_WallEnd;                   // for run details file
    clock_t                                            m_ClockStart;                // for run details file

    // member functions

    bool IsValidId(const int p_LogfileId)    { return ((p_LogfileId >= 0) && ((unsigned int)p_LogfileId < m_Logfiles.size())); }
    bool IsActiveId(const int p_LogfileId)   { return IsValidId(p_LogfileId) && m_Logfiles[p_LogfileId].active; }

    void ClearEntry(const int p_LogfileId) {
        if (IsValidId(p_LogfileId)) {
            m_Logfiles[p_LogfileId].active          = false;                        // not active
            m_Logfiles[p_LogfileId].logfiletype     = LOGFILE::NONE;
            m_Logfiles[p_LogfileId].filetype        = LOGFILETYPE::NONE;
            m_Logfiles[p_LogfileId].name            = "";
            m_Logfiles[p_LogfileId].timestamp       = false;
            m_Logfiles[p_LogfileId].label           = false;
            m_Logfiles[p_LogfileId].h5File.fileId   = -1;
            m_Logfiles[p_LogfileId].h5File.groupId  = -1;
            m_Logfiles[p_LogfileId].h5File.dataSets = {};
        }
    }

    bool DoIt(const string p_Class, const int p_Level, const std::vector<string> p_EnabledClasses, const int p_EnabledLevel);
    void Say_(const string p_SayStr);
    bool Write_(const int p_LogfileId, const string p_LogStr);
    bool Write_(const int p_LogfileId, const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues, const bool p_Flush = false);
    bool WriteHDF5_(h5AttrT& p_H5file, const string p_H5filename, const size_t p_DataSetIdx);
    bool Flush_(const int p_LogfileId) { return Write_(p_LogfileId, {}, true); }
    bool Put_(const int p_LogfileId, const string p_LogStr, const string p_Label = "");
    bool Put_(const int p_LogfileId, const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues);
    bool Debug_(const string p_DbgStr);
    bool Close_(const int p_LogfileId);

    PROPERTY_DETAILS StellarPropertyDetails(const ANY_STAR_PROPERTY p_Property);
    PROPERTY_DETAILS BinaryPropertyDetails(const BINARY_PROPERTY p_Property);
    PROPERTY_DETAILS ProgramOptionDetails(const PROGRAM_OPTION p_Property, const size_t p_Idx = 0);
    STR_STR_STR_STR  FormatFieldHeaders(const PROPERTY_DETAILS p_Details, string p_HeaderSuffix = "");
    LogfileDetailsT  StandardLogFileDetails(const LOGFILE p_Logfile, const string p_FileSuffix = "");

    std::tuple<bool, LOGFILE> GetLogfileDescriptorKey(const string p_Value);
    std::tuple<bool, LOGFILE> GetStandardLogfileKey(const int p_FileId);

    bool  OpenHDF5RunDetailsFile(const string p_Filename = RUN_DETAILS_FILE_NAME);
    hid_t CreateHDF5Dataset(const string p_Filename, const hid_t p_GroupId, const string p_DatasetName, const hid_t p_H5DataType, const string p_UnitsStr, const size_t p_HDF5ChunkSize);
    hid_t GetHDF5DataType(const TYPENAME p_COMPASdatatype, const int p_FieldWidth, const STRING_QUALIFIER p_StringQualifier = STRING_QUALIFIER::FIXED_LENGTH);

    bool NotesPropertyPresent(const ANY_PROPERTY_VECTOR p_RecordProperties) { return std::find(p_RecordProperties.begin(), p_RecordProperties.end(), T_ANY_PROPERTY(PROGRAM_OPTION::NOTES)) != p_RecordProperties.end(); }


    /*
     * Construct a log record to be written to a standard logfile
     *
     * This function constructs a log record to be written to one of the standard COMPAS logfiles.
     * The record to be constructed is identified by the logfile to which it is to be written, and the data
     * is assembled on-the-fly - except possibly for the parameter p_SpecifiedProperty (see below).  The 
     * star from which the data should be gathered is passed as a parameter, as is the logfile for which the
     * record should be constructed.  The logfile record properties and format vector appropriate to the
     * specified logfile are also passed as parameters - this function does not need to retrieve/construct
     * them.
     *
     * This function will, if the value of parameter p_UseSpecifiedValue is true, replace the value of any 
     * property in the record matching the value of the parameter p_SpecifiedProperty with the value specified
     * by the parameter p_SpecifiedPropertyValue.  (See the variant function GetLogStandardRecord() for an
     * explanation of why p_UseSpecifiedValue is required).
     * 
     *
     * template <class T1, typename T2>
     * std::tuple<std::string, std::vector<COMPAS_VARIABLE_TYPE>> GetLogStandardRecord(const LOGFILE             p_LogFile,
     *                                                                                 const LOGRECORDTYPE       p_RecordType,
     *                                                                                 const T1* const           p_Star,
     *                                                                                 const ANY_PROPERTY_VECTOR p_RecordProperties,
     *                                                                                 const std::vector<string> p_FmtVector,
     *                                                                                 const std::vector<bool>   p_Annotations,
     *                                                                                 const bool                p_UseSpecifiedValue,
     *                                                                                 const ANY_STAR_PROPERTY   p_SpecifiedProperty,
     *                                                                                 const T2                  p_SpecifiedPropertyValue)
     *
     * @param   [IN]    p_LogFile                   The logfile for which the record should be constructed
     * @param   [IN]    p_RecordType                The logfile record type
     * @param   [IN]    p_Star                      The object from which the values to be used to construct the record should be retrieved
     * @param   [IN]    p_RecordProperties          The logfile record properties pertaining to p_LogFile
     * @param   [IN]    p_FmtVector                 The logfile format vector pertaining to p_LogFile
     * @param   [IN]    p_Annotations               The logfile annotations vector pertaining to p_LogFile
     * @param   [IN]    p_UseSpecifiedValue         Flag to specify whether the parameters p_SpecifiedProperty and p_SpecifiedPropertyValue
     *                                              should be used to replace values retrieved from p_Star
     * @param   [IN]    p_SpecifiedProperty         The property type of the value to be replaced by p_SpecifiedPropertyValue
     * @param   [IN]    p_SpecifiedPropertyValue    The value of the property to be replaced
     * 
     * @return                                      tuple containing
     *                                                  - String formatted as log file record - empty string if an error occurred
     *                                                  - Vector of property values - empty vector if an error occurred
     */
    template <class T1, typename T2>
    std::tuple<std::string, std::vector<COMPAS_VARIABLE_TYPE>> GetLogStandardRecord(const LOGFILE             p_LogFile,
                                                                                    const LOGRECORDTYPE       p_RecordType,
                                                                                    const T1* const           p_Star,
                                                                                    const ANY_PROPERTY_VECTOR p_RecordProperties,
                                                                                    const std::vector<string> p_FmtVector,
                                                                                    const std::vector<bool>   p_Annotations,
                                                                                    const bool                p_UseSpecifiedValue,
                                                                                    const ANY_STAR_PROPERTY   p_SpecifiedProperty,
                                                                                    const T2                  p_SpecifiedPropertyValue) {

        bool ok = true;                                                                                                         // initially

        bool hdf5 = OPTIONS->LogfileType() == LOGFILETYPE::HDF5;                                                                // logging to hdf5 file?

        // construct log record from current data

        string logRecord = "";                                                                                                  // for CSV, TSV, TXT files: the record to be written to the log file
        std::vector<COMPAS_VARIABLE_TYPE> logRecordValues = {};                                                                 // for HDF5 files: vector of values to be written
                                                             
        // set delimiter based on logfile type
        string delimiter = "";                                                                                                  // default
        switch (OPTIONS->LogfileType()) {
            case LOGFILETYPE::HDF5: delimiter = ""; break;                                                                      // HDF5
            case LOGFILETYPE::CSV : delimiter = DELIMITERValue.at(DELIMITER::COMMA); break;                                     // CSV
            case LOGFILETYPE::TSV : delimiter = DELIMITERValue.at(DELIMITER::TAB); break;                                       // TSV
            case LOGFILETYPE::TXT : delimiter = DELIMITERValue.at(DELIMITER::SPACE); break;                                     // TXT
            default               : delimiter = ""; break;                                                                      // default
        }

        // get values
        //    - format for printing for CSV, TSV and TXT files
        //    - record for HDF5 files
        COMPAS_VARIABLE_TYPE value;                                                                                             // property value
        string               valueStr;                                                                                          // string for formatted value

        int index = 0;
        for (auto &property : p_RecordProperties) {                                                                             // for each property to be included in the log record
            ANY_PROPERTY_TYPE thisPropertyType = boost::apply_visitor(VariantPropertyType(), property);                         // get property type for this property
            
            boost::variant<string> fmtStr(p_FmtVector[index++]);                                                                // format string for this property

            // program option NOTES is special...
            //
            // NOTES is stored as a vector of strings, and we want to print all notes specified by the
            // annotations specification for the file

            bool needPropertyValue = true;
            if (thisPropertyType == ANY_PROPERTY_TYPE::T_PROGRAM_OPTION ) {                                                     // this property a PROGRAM_OPTION property?
                PROGRAM_OPTION thisProperty = boost::get<PROGRAM_OPTION>(property);                                             // get property
                if (thisProperty == PROGRAM_OPTION::NOTES) {                                                                    // PROGRAM_OPTION::NOTES?
                                                                                                                                // yes
                    for (size_t idx = 0; idx < p_Annotations.size(); idx ++) {                                                  // for each user-specified annotation
                        if (p_Annotations[idx]) {                                                                               // include it?
                            value = boost::variant<string>(OPTIONS->Notes(idx));                                                // yes - get value

                            if (hdf5) {                                                                                         // HDF5 file?
                                logRecordValues.push_back(value);                                                               // yes - add value to vector of values
                            }
                            else {                                                                                              // no - CSV, TSV, or TXT file
                                valueStr = boost::apply_visitor(FormatVariantValue(), value, fmtStr);                           // format value
                                logRecord += valueStr + delimiter;                                                              // add value string to log record - with delimiter
                            }
                        }
                    }
                    needPropertyValue = false;                                                                                  // have property value (or don't need it)
                }
            }

            if (needPropertyValue) {                                                                                            // have property value yet?
                                                                                                                                // no - need to get it
                ANY_STAR_PROPERTY thisProperty = ANY_STAR_PROPERTY::NONE;
                if (p_UseSpecifiedValue) {                                                                                      // replace specific value?

                    // determine if this property is p_SpecifiedProperty - if it is, replace the value with p_SpecifiedPropertyValue
                    // (this is specific to replacing a STELLAR_TYPE property - may need to generalise this to other types later - cross that bridge then)

                    switch (thisPropertyType) {
                        case ANY_PROPERTY_TYPE::T_STAR_PROPERTY     : { STAR_PROPERTY      prop = boost::get<STAR_PROPERTY>(property);      thisProperty = (ANY_STAR_PROPERTY)prop; } break;
                        case ANY_PROPERTY_TYPE::T_STAR_1_PROPERTY   : { STAR_1_PROPERTY    prop = boost::get<STAR_1_PROPERTY>(property);    thisProperty = (ANY_STAR_PROPERTY)prop; } break;
                        case ANY_PROPERTY_TYPE::T_STAR_2_PROPERTY   : { STAR_2_PROPERTY    prop = boost::get<STAR_2_PROPERTY>(property);    thisProperty = (ANY_STAR_PROPERTY)prop; } break;
                        case ANY_PROPERTY_TYPE::T_SUPERNOVA_PROPERTY: { SUPERNOVA_PROPERTY prop = boost::get<SUPERNOVA_PROPERTY>(property); thisProperty = (ANY_STAR_PROPERTY)prop; } break;
                        case ANY_PROPERTY_TYPE::T_COMPANION_PROPERTY: { COMPANION_PROPERTY prop = boost::get<COMPANION_PROPERTY>(property); thisProperty = (ANY_STAR_PROPERTY)prop; } break;
                        default: ok = false;                                                                                    // unknown property type - that's not ok...
                    }
                }

                // if not ok at this point we'll just take the current value - if there is really a problem
                // that will fail (as it would have normally if we weren't trying to replace a value)

                if (ok && p_UseSpecifiedValue && (thisProperty == p_SpecifiedProperty)) {                                       // replace specified property?
                    value    = p_SpecifiedPropertyValue;                                                                        // yes - use value passed as parameter
                    if (hdf5) {                                                                                                 // yes - HDF5 file?
                        logRecordValues.push_back(value);                                                                       // yes - add value to vector of values
                    }
                    else {                                                                                                      // no - CSV, TSV, or TXT file
                        valueStr = boost::apply_visitor(FormatVariantValue(), value, fmtStr);                                   // format value
                        logRecord += valueStr + delimiter;                                                                      // add value string to log record - with delimiter
                    }
                }
                else {                                                                                                          // use current value
                    std::tie(ok, value) = p_Star->PropertyValue(property);                                                      // get property flag and value
                    if (ok) {                                                                                                   // have valid property value
                        if (hdf5) {                                                                                             // yes - HDF5 file?
                            logRecordValues.push_back(value);                                                                   // yes - add value to vector of values
                        }
                        else {                                                                                                  // no - CSV, TSV, or TXT file
                            valueStr = boost::apply_visitor(FormatVariantValue(), value, fmtStr);                               // format value
                            logRecord += valueStr + delimiter;                                                                  // add value string to log record - with delimiter
                        }
                    }
                    else {                                                                                                      // unknown property type - should never happen
                        Squawk(ERR_MSG(ERROR::UNKNOWN_PROPERTY_TYPE));                                                          // show warning
                        logRecord = "";                                                                                         // empty record
                        logRecordValues.clear();                                                                                // and values vector
                        break;                                                                                                  // stop now
                    }
                }
                ok = true;
            }
        }

        if (ok) {

            // if we are writing to the SSE Switch file we add two pre-defined columns to 
            // the log record.  These are:
            //
            // ( i) the stellar type from which the star is switching
            // (ii) the stellar type to which the star is switching
            //
            // if we are writing to the BSE Switch file we add three pre-defined columns
            // to the end of the log record.  These are:
            //
            // (  i) the star switching - 1 = primary, 2 = secondary
            // ( ii) the stellar type from which the star is switching
            // (iii) the stellar type to which the star is switching
            //
            // These are hard-coded here rather than in the *_PROPERTY_DETAIL maps in
            // constants.h so that they will always be present in the switch file -
            // this way users can't add or remove them at runtime via the logfile-definitions
            // option.

            string fmtStr = "%4.1d";                                                                                            // format - all integers here

            if (p_LogFile == LOGFILE::BSE_SWITCH_LOG) {
                int starSwitching = m_PrimarySwitching ? 1 : 2;                                                                 // primary (1) or secondary (2)
                if (hdf5) {                                                                                                     // yes - HDF5 file?
                    logRecordValues.push_back(starSwitching);                                                                   // add value to vector of values
                }
                else {                                                                                                          // no - CSV, TSV, or TXT file
                    logRecord += utils::vFormat(fmtStr.c_str(), starSwitching) + delimiter;                                     // add value string to log record - with delimiter
                }
            }

            if (p_LogFile == LOGFILE::BSE_SWITCH_LOG || p_LogFile == LOGFILE::SSE_SWITCH_LOG) {
                if (hdf5) {                                                                                                     // HDF5 file?
                    logRecordValues.push_back(m_TypeSwitchingFrom);                                                             // yes - add value to vector of values
                }
                else {                                                                                                          // no - CSV, TSV, or TXT file
                    logRecord += utils::vFormat(fmtStr.c_str(), m_TypeSwitchingFrom) + delimiter;                               // add value string to log record - with delimiter
                }

                if (hdf5) {                                                                                                     // HDF5 file?
                    logRecordValues.push_back(m_TypeSwitchingTo);                                                               // yes - add value to vector of values
                }
                else {                                                                                                          // no - CSV, TSV, or TXT file
                    logRecord += utils::vFormat(fmtStr.c_str(), m_TypeSwitchingTo) + delimiter;                                 // add value string to log record - with delimiter
                }
            }

            // we add the record type column to the end of the log record here for all logfiles
            // except the switch files (BSE_SWITCH_LOG and SSE_SWOTCH_LOG).
            //
            // This is hard-coded here rather than in the *_PROPERTY_DETAIL maps in constants.h
            // so that it will always be present in the logfile - this way users can't add or 
            // remove it at runtime via the logfile-definitions option.

            if (p_LogFile != LOGFILE::BSE_SWITCH_LOG && p_LogFile != LOGFILE::SSE_SWITCH_LOG) {                                 // switch file?
                fmtStr = "%10.1u";                                                                                              // no - proceed
                if (hdf5) {                                                                                                     // HDF5 file?
                    logRecordValues.push_back(p_RecordType);                                                                    // add value to vector of values
                }
                else {                                                                                                          // no - CSV, TSV, or TXT file
                    logRecord += utils::vFormat(fmtStr.c_str(), p_RecordType);                                                  // add value string to log record
                }
            }
        }

        return std::make_tuple(logRecord, logRecordValues);
    }


    /*
     * This variant of GetLogStandardRecord() is here because I can't readily figure out how to 
     * have an optional parameter in a template function if not providing that optional parameter 
     * doesn't give the compiler enough context to deduce the type of the template type.
     * 
     * Clear as mud?  Well, for the declaration of GetLogStandardRecord() above, template type T2
     * can't be deduced if the parameter p_SpecifiedPropertyValue is not supplied - even if it
     * is given a default value in the parameter list.  There must be a way of making it work, but
     * I can't see it - and I don't want to spend any more time trying to figure it out.  Maybe 
     * later.  In the meantime, this method works.
     * 
     * See description of GetLogStandardRecord() for functionality and parameter descriptions.
     * 
     */
    template <class T>
    std::tuple<string, std::vector<COMPAS_VARIABLE_TYPE>> GetLogStandardRecord(const LOGFILE             p_LogFile,
                                                                               const LOGRECORDTYPE       p_RecordType,
                                                                               const T* const            p_Star,
                                                                               const ANY_PROPERTY_VECTOR p_RecordProperties,
                                                                               const std::vector<string> p_FmtVector,
                                                                               const std::vector<bool>   p_Annotations) {

        return GetLogStandardRecord(p_LogFile, 
                                    p_RecordType,
                                    p_Star, 
                                    p_RecordProperties, 
                                    p_FmtVector, 
                                    p_Annotations,
                                    false, 
                                    ANY_STAR_PROPERTY::STELLAR_TYPE, 
                                    STELLAR_TYPE::NONE);
    }


    /*
     * Write a record to a standard logfile
     *
     * This function writes a log record to one of the standard COMPAS logfiles.
     * The record to be written is identified by the logfile to which it is to be written, and the data
     * is assembled on-the-fly.  The star from which the data should be gathered is passed as a parameter,
     * as is the logfile to which the record should be written, as well as a string suffix for the logfile
     * name.
     *
     * If the logfile is not already open, it will be opened by calling this function.
     *
     * The logging classes and level set by program options are honoured by this function.  The log class
     * and level for this record are passed as parameters.  If the log class for this record is not
     * configured to be logged, this function returns without performing any logging action (i.e. the log
     * file status (open, closed) and contents, will be unchanged by this function.  Similarly, if the
     * log level for this record indicates that the record should not be logged, no changes will be made
     * by this function.
     * 
     *
     * template <class T>
     * bool LogStandardRecord(const string        p_LogClass,
     *                        const int           p_LogLevel,
     *                        const LOGFILE       p_LogFile,
     *                        const LOGRECORDTYPE p_RecordType,
     *                        const T* const      p_Star,
     *                        const string        p_FileSuffix = "")
     *
     * @param   [IN]    p_LogClass                  Class to determine if record should be written
     * @param   [IN]    p_RecordType                The logfile record type
     * @param   [IN]    p_LogLevel                  Level to determine if record should be written
     * @param   [IN]    p_LogFile                   The logfile to which the record should be written
     * @param   [IN]    p_Star                      The star object from which the field values should be retrieved
     * @param   [IN]    p_FileSuffix                String suffix to be added to the logfile name (optional, default = "")
     * @return                                      Boolean status (true = success, false = failure)
     */
    template <class T>
    bool LogStandardRecord(const string        p_LogClass,
                           const int           p_LogLevel,
                           const LOGFILE       p_LogFile,
                           const LOGRECORDTYPE p_RecordType,
                           const T* const      p_Star,
                           const string        p_FileSuffix = "") {

        bool ok = true;                                                                                                     // initially

        LogfileDetailsT fileDetails = StandardLogFileDetails(p_LogFile, p_FileSuffix);                                      // get record details - open file (if necessary)
        if (fileDetails.id >= 0) {                                                                                          // file open?
            if (((1 << (p_RecordType - 1)) & fileDetails.recordTypes) > 0) {                                                // yes - record type enabled?
                                                                                                                            // yes - proceed
                string logRecordString;                                                                                     // for CSV, TSV, TXT files: the record to be written to the log file
                std::vector<COMPAS_VARIABLE_TYPE> logRecordValues;                                                          // for HDF5 files: vector of values to be written

                // construct the record - gets both string and vector of values
                std::tie(logRecordString, logRecordValues) = GetLogStandardRecord(p_LogFile, p_RecordType, p_Star, fileDetails.recordProperties, fileDetails.fmtStrings, fileDetails.annotations);

                if (OPTIONS->LogfileType() == LOGFILETYPE::HDF5)                                                            // logging to HDF5 file?
                    ok = Put(fileDetails.id, p_LogClass, p_LogLevel, logRecordValues);                                      // yes - write the record
                else                                                                                                        // not HDF5
                    ok = Put(fileDetails.id, p_LogClass, p_LogLevel, logRecordString);                                      // write the record

                if (!ok) Squawk(ERR_MSG(ERROR::FILE_WRITE_ERROR) + " while writing to logfile " + fileDetails.filename);    // show warning if record not written ok
            }
        }
        return ok;
    }


    /*
     * The following two variants of LogStandardRecord() are here to allow logging records already
     * constructed.  Parameters are the same except for an additional parameter in each case: for
     * HDF5 files the extra parameter is a vector of values, and for non-HDF5 files the parameter
     * is a string.
     * 
     * See description of GetLogStandardRecord() for functionality and parameter descriptions.
     */

    /*
     * ... for non-HDF5 files
     * @param   [IN]    p_LogRecordString           The previously constructed string to be written to the file
     */
    template <class T>
    bool LogStandardRecord(const string        p_LogClass,
                           const int           p_LogLevel,
                           const LOGFILE       p_LogFile,
                           const LOGRECORDTYPE p_RecordType,
                           const T* const      p_Star,
                           const string        p_FileSuffix,
                           const string        p_LogRecordString) {

        bool ok = true;                                                                                                     // initially

        LogfileDetailsT fileDetails = StandardLogFileDetails(p_LogFile, p_FileSuffix);                                      // get record details - open file (if necessary)
        if (fileDetails.id >= 0) {                                                                                          // file open?
            if ((p_RecordType & fileDetails.recordTypes) > 0) {                                                             // yes - record type enabled?
                                                                                                                            // yes - proceed
                if (m_Logfiles[fileDetails.id].filetype == LOGFILETYPE::HDF5) {                                             // HDF5 file?
                    Squawk(ERR_MSG(ERROR::UNEXPECTED_LOG_FILE_TYPE) + " while writing to logfile " + fileDetails.filename); // yes - show warning: unexpected logfile type
                    ok = false;                                                                                             // fail
                }
                else {                                                                                                      // not HDF5
                    ok = Put(fileDetails.id, p_LogClass, p_LogLevel, p_LogRecordString);                                    // write the record
                    if (!ok) Squawk(ERR_MSG(ERROR::FILE_WRITE_ERROR) + " while writing to logfile " + fileDetails.filename);// show warning if record not written ok
                }
            }
        }
        return ok;
    }

    /*
     * ... for HDF5 files
     * @param   [IN]    p_LogRecordValues           The previously constructed string to be written to the file
     */
    template <class T>
    bool LogStandardRecord(const string                            p_LogClass,
                           const int                               p_LogLevel,
                           const LOGFILE                           p_LogFile,
                           const LOGRECORDTYPE                     p_RecordType,
                           const T* const                          p_Star,
                           const string                            p_FileSuffix,
                           const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues) {

        bool ok = true;                                                                                                     // initially

        LogfileDetailsT fileDetails = StandardLogFileDetails(p_LogFile, p_FileSuffix);                                      // get record details - open file (if necessary)
        if (fileDetails.id >= 0) {                                                                                          // file open?
            if ((p_RecordType & fileDetails.recordTypes) > 0) {                                                             // yes - record type enabled?
                                                                                                                            // yes - proceed
                if (m_Logfiles[fileDetails.id].filetype == LOGFILETYPE::HDF5) {                                             // HDF5 file?
                    ok = Put(fileDetails.id, p_LogClass, p_LogLevel, p_LogRecordValues);                                    // yes - write the record
                    if (!ok) Squawk(ERR_MSG(ERROR::FILE_WRITE_ERROR) + " while writing to logfile " + fileDetails.filename);// show warning if record not written ok
                }
                else {                                                                                                      // not HDF5
                    Squawk(ERR_MSG(ERROR::UNEXPECTED_LOG_FILE_TYPE) + " while writing to logfile " + fileDetails.filename); // show warning: unexpected logfile type
                    ok = false;                                                                                             // fail
                }
            }
        }
        return ok;
    }


    void PrintLogfileRecordDetails(const ANY_PROPERTY_VECTOR& p_LogfileRecord, const string p_LogfileRecordName);

    void UpdateLogfileRecordSpecs(const LOGFILE             p_Logfile,
                                  bool                      p_UseDefaultProps,
                                  const ANY_PROPERTY_VECTOR p_AddProps,
                                  const ANY_PROPERTY_VECTOR p_SubtractProps,
                                  const std::vector<bool>   p_AddNotes,
                                  const std::vector<bool>   p_SubtractNotes);

    bool UpdateAllLogfileRecordSpecs();


public:

    ~Log() {                                                                                                                            // destructor
        Stop();                                                                                                                         // stop logging - flushes and closes all open logfiles
        delete m_Instance;                                                                                                              // delete the instance variable
    }

    // instance function
    static Log* Instance();                                                                                                             // the singleton instance is exposed, not the constructor


    // member functions
    void   Start(const string              p_LogBasePath,
                 const string              p_LogContainerName,
                 const string              p_LogNamePrefix,
                 const int                 p_LogLevel,
                 const std::vector<string> p_LogClasses,
                 const int                 p_DbgLevel,
                 const std::vector<string> p_DbgClasses,
                 const bool                p_DbgToFile,
                 const bool                p_ErrorsToFile,
                 const LOGFILETYPE         p_LogfileType);

    void   Stop(std::tuple<int, int> p_ObjectStats = std::make_tuple(0, 0));

    bool   Enabled() const { return m_Enabled; }

    int    Open(const string p_LogFileName, const bool p_Append, const bool p_TimeStamp, const bool p_Label, const LOGFILE p_StandardLogfile = LOGFILE::NONE);
    bool   Close(const int p_LogfileId);

    bool   Write(const int p_LogfileId, const string p_LogClass, const int p_LogLevel, const string p_LogStr);
    bool   Write(const int p_LogfileId, const string p_LogClass, const int p_LogLevel, const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues, const bool p_Flush = false);
    
    bool   Put(const int p_LogfileId, const string p_LogClass, const int p_LogLevel, const string p_LogStr);
    bool   Put(const int p_LogfileId, const string p_LogClass, const int p_LogLevel, const std::vector<COMPAS_VARIABLE_TYPE> p_LogRecordValues);

    bool   Debug(const string p_DbgClass, const int p_DbgLevel, const string p_DbgStr);
    bool   DebugWait(const string p_DbgClass, const int p_DbgLevel, const string p_DbgStr);

    bool   Error(const string p_ErrStr);

    void   Squawk(const string squawkStr);

    void   Say(const string p_SayClass, const int p_SayLevel, const string p_SayStr);


    // SetSwitchParameters is called by Star::SwitchTo to set the parameters 
    // to be written to the BSE Switch Log file
    void   SetSwitchParameters(const OBJECT_ID    p_ObjectIdSwitching, 
                               const STELLAR_TYPE p_TypeSwitchingFrom, 
                               const STELLAR_TYPE p_TypeSwitchingTo) {
        m_ObjectIdSwitching = p_ObjectIdSwitching;                          // the object id of the Star object switching stellar type
        m_TypeSwitchingFrom = p_TypeSwitchingFrom;                          // the stellar type from which the Star object is switching
        m_TypeSwitchingTo   = p_TypeSwitchingTo;                            // the stellar type to which the Star object is switching
    }

    OBJECT_ID ObjectIdSwitching() { return m_ObjectIdSwitching; }



    // standard logfile logging functions

    bool CloseStandardFile(const LOGFILE p_LogFile, const bool p_Erase = true);
    bool CloseAllStandardFiles();

    std::tuple<ANY_PROPERTY_VECTOR, std::vector<string>, std::vector<bool>> GetStandardLogFileRecordDetails(const LOGFILE p_Logfile);

    template <class T>
    bool LogBeBinary(const T* const p_Binary,
                     const BE_BINARY_RECORD_TYPE p_RecordType)                      { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_BE_BINARIES)), 0, LOGFILE::BSE_BE_BINARIES, static_cast<LOGRECORDTYPE>(p_RecordType), p_Binary); }

    template <class T>
    bool LogBSEDetailedOutput(const T* const p_Binary, 
                              const long int p_Id,
                              const BSE_DETAILED_RECORD_TYPE p_RecordType)          { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_DETAILED_OUTPUT)), 0, LOGFILE::BSE_DETAILED_OUTPUT, static_cast<LOGRECORDTYPE>(p_RecordType), p_Binary, "_" + std::to_string(abs(p_Id))); }

    template <class T>
    bool LogBSEPulsarEvolutionParameters(const T* const p_Binary,
                                         const PULSAR_RECORD_TYPE p_RecordType)     { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_PULSAR_EVOLUTION)), 0, LOGFILE::BSE_PULSAR_EVOLUTION, static_cast<LOGRECORDTYPE>(p_RecordType), p_Binary); }

    template <class T>
    bool LogBSESupernovaDetails(const T* const p_Binary,
                                const BSE_SN_RECORD_TYPE p_RecordType)              { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SUPERNOVAE)), 0, LOGFILE::BSE_SUPERNOVAE, static_cast<LOGRECORDTYPE>(p_RecordType), p_Binary); }
    
    template <class T>
    bool LogBSESwitchLog(const T* const p_Binary, const bool p_PrimarySwitching) {
        m_PrimarySwitching = p_PrimarySwitching;        
        return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SWITCH_LOG)), 0, LOGFILE::BSE_SWITCH_LOG, 1U, p_Binary);
    }

    template <class T>
    bool LogBSESystemParameters(const T* const p_Binary,
                                const BSE_SYSPARMS_RECORD_TYPE p_RecordType)        { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_SYSTEM_PARAMETERS)), 0, LOGFILE::BSE_SYSTEM_PARAMETERS, static_cast<LOGRECORDTYPE>(p_RecordType), p_Binary); }

    template <class T>
    bool LogCommonEnvelope(const T* const p_Binary,
                           const CE_RECORD_TYPE p_RecordType)                       { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_COMMON_ENVELOPES)), 0, LOGFILE::BSE_COMMON_ENVELOPES, static_cast<LOGRECORDTYPE>(p_RecordType), p_Binary); }

    template <class T>
    bool LogDoubleCompactObject(const T* const p_Binary,
                                const DCO_RECORD_TYPE p_RecordType)                 { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS)), 0, LOGFILE::BSE_DOUBLE_COMPACT_OBJECTS, static_cast<LOGRECORDTYPE>(p_RecordType), p_Binary); }

    template <class T>
    bool LogRLOFParameters(const T* const p_Binary,
                           const RLOF_RECORD_TYPE p_RecordType)                     { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::BSE_RLOF_PARAMETERS)), 0, LOGFILE::BSE_RLOF_PARAMETERS, static_cast<LOGRECORDTYPE>(p_RecordType), p_Binary); }

    template <class T>
    bool LogSSEDetailedOutput(const T* const p_Star, const int p_Id,
                              const SSE_DETAILED_RECORD_TYPE p_RecordType)          { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_DETAILED_OUTPUT)), 0, LOGFILE::SSE_DETAILED_OUTPUT, static_cast<LOGRECORDTYPE>(p_RecordType), p_Star, "_" + std::to_string(abs(p_Id))); }

    template <class T>
    bool LogSSESupernovaDetails(const T* const p_Star,
                                const SSE_SN_RECORD_TYPE p_RecordType)              { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SUPERNOVAE)), 0, LOGFILE::SSE_SUPERNOVAE, static_cast<LOGRECORDTYPE>(p_RecordType), p_Star); }

    template <class T>
    bool LogSSESwitchLog(const T* const p_Star)                                     { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SWITCH_LOG)), 0, LOGFILE::SSE_SWITCH_LOG, 1U, p_Star); }

    template <class T>
    bool LogSSESystemParameters(const T* const p_Star,
                                const SSE_SYSPARMS_RECORD_TYPE p_RecordType)        { return LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SYSTEM_PARAMETERS)), 0, LOGFILE::SSE_SYSTEM_PARAMETERS, static_cast<LOGRECORDTYPE>(p_RecordType), p_Star); }


    template <class T>
    void StashSSESupernovaDetails(const T* const p_Star, const STELLAR_TYPE p_StellarType, const SSE_SN_RECORD_TYPE p_RecordType) {

        m_SSESupernovae_DelayedWrite.logRecordType = static_cast<LOGRECORDTYPE>(p_RecordType);

        // if we don't already have the SSE Supernova log record properties that we need, get them
        // this will only need to be done once per run, so not a big overhead
        if (m_SSESupernovae_DelayedWrite.logRecordProperties.empty() || 
            m_SSESupernovae_DelayedWrite.logRecordFmtVector.empty()  || 
            m_SSESupernovae_DelayedWrite.logFileAnnotations.empty()) {

            std::tie(m_SSESupernovae_DelayedWrite.logRecordProperties, 
                     m_SSESupernovae_DelayedWrite.logRecordFmtVector, 
                     m_SSESupernovae_DelayedWrite.logFileAnnotations) = LOGGING->GetStandardLogFileRecordDetails(LOGFILE::SSE_SUPERNOVAE);
        }

        // get a formatted record with current data
        // this will replace any existing stashed record - no queue here
        std::tie(m_SSESupernovae_DelayedWrite.logRecordString, 
                 m_SSESupernovae_DelayedWrite.logRecordValues) = GetLogStandardRecord(LOGFILE::SSE_SUPERNOVAE,
                                                                                      m_SSESupernovae_DelayedWrite.logRecordType, 
                                                                                      p_Star, 
                                                                                      m_SSESupernovae_DelayedWrite.logRecordProperties, 
                                                                                      m_SSESupernovae_DelayedWrite.logRecordFmtVector,
                                                                                      m_SSESupernovae_DelayedWrite.logFileAnnotations, 
                                                                                      true, 
                                                                                      (ANY_STAR_PROPERTY)STAR_PROPERTY::STELLAR_TYPE, 
                                                                                      p_StellarType);
    }

    template <class T>
    bool LogStashedSSESupernovaDetails(const T* const p_Star) { 
        bool result = true;

        // if the stashed SSE Supernova record is non-empty, print it, then clear it - otherwise do nothing

        if (OPTIONS->LogfileType() == LOGFILETYPE::HDF5) {                  // logging to HDF5 file?
            if (!m_SSESupernovae_DelayedWrite.logRecordValues.empty()) {    // yes - need to log?
                result = LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SUPERNOVAE)), 0, LOGFILE::SSE_SUPERNOVAE, m_SSESupernovae_DelayedWrite.logRecordType, p_Star, "", m_SSESupernovae_DelayedWrite.logRecordValues);
                m_SSESupernovae_DelayedWrite.logRecordValues = {};          // clear record
            }
        }
        else {                                                              // no - not HDF5
            if (!m_SSESupernovae_DelayedWrite.logRecordString.empty()) {    // need to log?
                result = LogStandardRecord(std::get<2>(LOGFILE_DESCRIPTOR.at(LOGFILE::SSE_SUPERNOVAE)), 0, LOGFILE::SSE_SUPERNOVAE, m_SSESupernovae_DelayedWrite.logRecordType, p_Star, "", m_SSESupernovae_DelayedWrite.logRecordString);
                m_SSESupernovae_DelayedWrite.logRecordString = "";          // clear record
            }
        }

        return result;
    }
};

#endif // __Log_h__
