/**
 * LoxcodeR
 */

#include <Rcpp.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <assert.h>
#include "edlib.h"
#include <regex>

using namespace Rcpp;
using namespace std;


void printAlignment(const char* query, const char* target,
                    const unsigned char* alignment, const int alignmentLength,
                    const int position, const EdlibAlignMode modeCode, int editDistance) {

  int tIdx = -1;
  int qIdx = -1;
  if (modeCode == EDLIB_MODE_HW) {
    tIdx = position;
    for (int i = 0; i < alignmentLength; i++) {
      if (alignment[i] != EDLIB_EDOP_INSERT)
        tIdx--;
    }
  }
  for (int start = 0; start < alignmentLength; start += 75) {
    // target
    printf("T: ");
    int startTIdx = -1;
    for (int j = start; j < start + 75 && j < alignmentLength; j++) {
      if (alignment[j] == EDLIB_EDOP_INSERT)
        printf("-");
      else
        printf("%c", target[++tIdx]);
      if (j == start)
        startTIdx = tIdx;
    }
    printf(" (%d - %d)\n", max(startTIdx, 0), tIdx);

    // match / mismatch
    printf("   ");
    for (int j = start; j < start + 75 && j < alignmentLength; j++) {
      printf(alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
    }
    printf("\n");

    // query
    printf("Q: ");
    int startQIdx = qIdx;
    for (int j = start; j < start + 75 && j < alignmentLength; j++) {
      if (alignment[j] == EDLIB_EDOP_DELETE)
        printf("-");
      else
        printf("%c", query[++qIdx]);
      if (j == start)
        startQIdx = qIdx;
    }
    printf(" (%d - %d)\n\n", max(startQIdx, 0), qIdx);
  }
}



template<typename T>
void pop_front(std::vector<T>& vec)
{
  assert(!vec.empty());
  vec.erase(vec.begin());
}

void readFASTA(ifstream &F , std::vector<string> & lines)
{
  getline(F,lines[0]);getline(F,lines[1]);getline(F,lines[2]);getline(F,lines[3]);
}

//barcode elements;
std::map<string, int> ele_R1={{"ACTCCGCA",1},{"TCCAGAATTTGTAT",2},{"ACATCCAC",3},{"AAAGGAATTTCTCC",4},{"ATTTCCTC",5},{"GCCCGAATTTTTTC",6},{"GCTACTGG",7},{"ATGAGAATTTATGG",8},{"AACTAGAA",9},{"TGCAGAATTTCCTC",10},{"CGACACTT",11},{"AACGGAATTTTCAA",12},{"CGTGTTTG",13},{"CAAACACG",-13},{"TTGAAAATTCCGTT",-12},{"AAGTGTCG",-11},{"GAGGAAATTCTGCA",-10},{"TTCTAGTT",-9},{"CCATAAATTCTCAT",-8},{"CCAGTAGC",-7},{"GAAAAAATTCGGGC",-6},{"GAGGAAAT",-5},{"GGAGAAATTCCTTT",-4},{"GTGGATGT",-3},{"ATACAAATTCTGGA",-2},{"TGCGGAGT",-1}};
std::map<string, string> ele_R2={{"GGATGAATTCGTGT","end"},{"ACTCCGCA","-1"},{"TCCAGAATTTGTAT","-2"},{"ACATCCAC","-3"},{"AAAGGAATTTCTCC","-4"},{"ATTTCCTC","-5"},{"GCCCGAATTTTTTC","-6"},{"GCTACTGG","-7"},{"ATGAGAATTTATGG","-8"},{"AACTAGAA","-9"},{"TGCAGAATTTCCTC","-10"},{"CGACACTT","-11"},{"AACGGAATTTTCAA","-12"},{"CGTGTTTG","-13"},{"CAAACACG","13"},{"TTGAAAATTCCGTT","12"},{"AAGTGTCG","11"},{"GAGGAAATTCTGCA","10"},{"TTCTAGTT","9"},{"CCATAAATTCTCAT","8"},{"CCAGTAGC","7"},{"GAAAAAATTCGGGC","6"},{"GAGGAAAT","5"},{"GGAGAAATTCCTTT","4"},{"GTGGATGT","3"},{"ATACAAATTCTGGA","2"},{"TGCGGAGT","1"},{"GTGCAAATTCGAGC","start"}};

std::vector<int> Consensus(vector<int> a, vector<int> b)
  /**
   * Used if R1, R2 overlap
   */
{
  int offset=0; bool do_align=false;
  int max_overlap=0,max_offset=0;
  int asize=a.size(),bsize=b.size();

  for(offset=-bsize; offset<asize; offset++) //check all possible starting positions of b
  {
    int match=0,not_match=0;

    for(int k=0; k<b.size();k++)
    {
      if(k+offset<0) ;
      else if(k+offset>=asize) ;
      else
      {
        if(a[k+offset]==b[k] && b[k]!=0) match++;
        if(a[k+offset]!=b[k] && a[k+offset]!=0 && b[k]!=0) not_match++;
      }
    }
    if(match>max_overlap && not_match==0) {max_overlap=match; max_offset=offset;}
  }

  std::vector<int> consensus,empty;
  if(max_overlap==0)
  {
    //add special case of full length
    if(a.size()+b.size()==12 && count(a.begin(),a.end(),0)==0 && count(b.begin(),b.end(),0)==0)
    {
      vector<int> full_length(a);
      full_length.push_back(0);
      full_length.insert(full_length.end(), b.begin(), b.end());
      //cout<<"full length"<<endl;
      return full_length;
    }
    else return empty;
  }


  for(int i=0; i<max_offset; i++) consensus.push_back(a[i]);
  for(int k=0; k<b.size();k++)
  {
    if(k+max_offset<a.size())
    {
      if(a[k+max_offset]==b[k]) {consensus.push_back(b[k]); continue;}
      if(a[k+max_offset]!=0 && b[k]==0){consensus.push_back(a[k+max_offset]); continue;}
      if(a[k+max_offset]==0 && b[k]!=0){consensus.push_back(b[k]); continue;}
    }
    else consensus.push_back(b[k]);
  }

  if(max_offset<0) { for(int i=max_offset;i<0;i++) pop_front(consensus); }


  return consensus;
}
//////////////////////////////
int H(string s1, string s2)
{
  int sum=0; for(int i=0 ;i<s1.size();i++) if(s1[i]!=s2[i]) sum++;
  return sum;
}
//////////////////////////////
int ED(string s1, string s2)
{
  return (edlibAlign(&s1[0u], s1.size(),&s2[0u], s2.size(), edlibNewAlignConfig(10, EDLIB_MODE_HW,  EDLIB_TASK_DISTANCE, NULL, 0)).editDistance);
}
//////////////////////////////
vector<int> map_lox(string R1, int start_loc)
{
  map<int,int> locs={{1,8},{43,14},{91,8},{133,14},{181,8},{223,14},{271,8}};
  vector<int> elements;
  int counter=0;
  for(auto l : locs){

    if(l.first+start_loc<0 || l.first+start_loc>R1.length()) {break;}
    string el=R1.substr(l.first+start_loc,l.second);
    auto it = ele_R1.find(el); // exact
    if(it!=ele_R1.end() && it->first.size() == l.second) {elements.push_back(it->second);  continue;}
    else
    {
      int inf_el=0; int min_h=1000, rep=0;
      for(auto i : ele_R1)
      {
        if(i.first.size()!=el.size()) continue;
        int h=H(i.first,el);
        if(h<=min_h){if(h==min_h) rep++; min_h=h; inf_el=i.second;}
      };

      if(min_h<=2 && rep==1) {elements.push_back(inf_el); cout<<"hamming "<<min_h<<endl;}
      else elements.push_back(0);
    }

    counter++;
  }

  return elements;
}

/////////////////////
int findNeedle(string A, string B)
{
  EdlibAlignResult result = edlibAlign(&A[0u],A.size(),&B[0u],B.size(), edlibNewAlignConfig(3, EDLIB_MODE_HW,  EDLIB_TASK_PATH, NULL, 0));
  int d=result.editDistance;
  //if(result.alignmentLength>0)printAlignment(&A[0u], &B[0u], result.alignment, result.alignmentLength, result.endLocations[0], EDLIB_MODE_HW, result.editDistance);
  edlibFreeAlignResult(result);
  return d;
}
/////////////////////
////////////////////
int findNeedlePos(string A, string B)
{
  EdlibAlignResult result = edlibAlign(&A[0u],A.size(),&B[0u],B.size(), edlibNewAlignConfig(3, EDLIB_MODE_HW,  EDLIB_TASK_PATH, NULL, 0));
  int d=result.startLocations[0];
  //if(d==5) printAlignment(&A[0u], &B[0u], result.alignment, result.alignmentLength, result.endLocations[0], EDLIB_MODE_HW, result.editDistance);
  edlibFreeAlignResult(result);
  return d;
}

//////////////////////
//////////////////////
vector<int> decR(string R, string head_seq, int head_length, string tail_seq, int tail_length, int i, int &all_mapped, int &ok_reads)
{
  vector<int> elements;
  int reads_missing_start=0;
  int tot_reads=i;

  string HEAD=R.substr(0,head_length);
  string TAIL=R.substr(R.size()-tail_length+1,tail_length);

  string loxF="ATAACTTCGTATAATGTATGCTATACGAAGTTAT",loxR="ATAACTTCGTATAGCATACATTATACGAAGTTAT";

  string start=head_seq; //"TCTAGAGGATCCCCGGGTACCGAGCTCGAATTTGCACATAACTTCGTATAATGTATGCTATACGAAGTTAT";
  string end=tail_seq; //"ATAACTTCGTATAGCATACATTATACGAAGTTATACACGAATTCATCCAGGC";

  int start_loc=0,end_loc=0;
  ///using edit ditance to find start (in read 1) and end (in read 2)
  // and we allow up to a single mismatch
  EdlibAlignResult  result1 = edlibAlign(&start[0u],start.size(),&HEAD[0u], HEAD.size(), edlibNewAlignConfig(10, EDLIB_MODE_HW,  EDLIB_TASK_PATH, NULL, 0));

  if(result1.alignmentLength>0)//first check if it aligns to full head
  {
    // printAlignment(&start[0u], &HEAD[0u], result1.alignment, result1.alignmentLength, result1.endLocations[0], EDLIB_MODE_HW, result1.editDistance);
    end_loc=result1.endLocations[0];

  }

  else //if not we try to rescue by aligning to LoxcodeF
  {
    EdlibAlignResult result2 = edlibAlign(&loxF[0u],loxF.size(),&HEAD[0u],HEAD.size(), edlibNewAlignConfig(5, EDLIB_MODE_HW,  EDLIB_TASK_PATH, NULL, 0));
    if(result2.alignmentLength>0)
    {
      //printAlignment(&loxF[0u], &HEAD[0u], result2.alignment, result2.alignmentLength, result2.endLocations[0], EDLIB_MODE_HW, result2.editDistance);
      end_loc=result2.endLocations[0];
    }
    edlibFreeAlignResult(result2);
  }

  edlibFreeAlignResult(result1);

  if(end_loc == 0){reads_missing_start++; cout<<"no start"<<endl; }

  else
  {
    ok_reads++;
    int max_r=R.size();
    if(findNeedle(end, TAIL)!=-1) {max_r=R.size()-tail_length + findNeedlePos(end, TAIL)+10; cout<<"found_end "<<max_r<<endl;}
    elements= map_lox(R.substr(0,max_r), end_loc);
  }
  /////////////////////////////////////////////
  return elements;
}


//' Decode FASTQ
//'
//' Recover loxcodes from raw Illumina FASTQ output
//' @param r Paths of R1, 2 respectively
//' @param meta User-defined data-frame for sample metadata
//' @param min_read_length min read length for R1, R2 filter respectively
//' @param full whether to supply full output (including read IDs, etc)
//' @param saturation whether to keep saturation information
//' @return S4 loxcode_sample object with decoded results
//' @export
// [[Rcpp::export]]
Rcpp::S4 decode(std::vector<std::string> r, std::string name, Rcpp::DataFrame meta,
                int min_r1_len, int min_r2_len, bool full, bool sat){

  ifstream fileR1(r[0]); ifstream fileR2(r[1]); // input files
  int counter=0;
  /*
   * code_readout.first = cassette sequence
   * code_readout.second = line numbers in fastq files corresponding to each
   */
  map< std::vector<int>, std::vector<int> > code_readout;
  vector<int> saturation; // track saturation

  int tot_reads = 0;
  // breakdown of discarded reads
  int too_small_reads = 0;
  int reads_missing_start = 0;
  int reads_missing_end = 0;
  int reads_multi_start = 0;
  int reads_multi_end = 0;
  int reads_consensus_filtered = 0;

  std::vector<std::vector<std::string> > reads_consensus_filtered_data;

  int keep=0; //keep track of reads that are useful

  if(!fileR1.is_open() || !fileR2.is_open()){
    // file not found
    Rcpp::stop("no such file");
  }

  double usable_reads=0;
  int ok_reads=0,all_mapped=0;

  for(int i=0; ;i++){
    // i is the read counter (starts from zero)
    if(fileR1.eof() || fileR2.eof()) break;

    std::vector<string> lines(4);
    readFASTA(fileR1, lines);
    string R1=lines[1];
    readFASTA(fileR2, lines);
    string R2=lines[1];

    if(R1.size()==0 || R2.size()==0) continue;

    tot_reads++;

    /////////////////////////////////////////////////

    if(R1.size()<75 || R2.size()<70) {cout<<"read too small"<<endl; continue;}
    cout<<endl<<"read # "<<tot_reads<<endl;

    cout<<"R1 ";
    string start_R1="TCTAGAGGATCCCCGGGTACCGAGCTCGAATTTGCACATAACTTCGTATAATGTATGCTATACGAAGTTAT";
    string end_R1="ATAACTTCGTATAGCATACATTATACGAAGTTATACACGAATTCATCCAGGC";
    vector<int> r1= decR(R1, start_R1, 75, end_R1, 60, i, all_mapped, ok_reads);
    for(auto k: r1) cout<<k<<" "; cout<<endl;

    string start_R2="GCCTGGATGAATTCGTGTATAACTTCGTATAATGTATGCTATACGAAGTTAT";
    string end_R2="GTGCAAATTCGAGCT";

    cout<<"R2 ";
    vector<int> r2 = decR(R2, start_R2, 60, end_R2, 70, i, all_mapped, ok_reads);
    std::reverse(r2.begin(),r2.end()); for(int j=0; j<r2.size(); j++) r2[j]*=-1;
    for(auto k: r2) cout<<k<<" "; cout<<endl;

    if(r1.size()!=r2.size())
    {
      cout<<"unequal"<<endl;
      for(auto k: r1) cout<<k<<" "; cout<<endl;
      for(auto k: r2) cout<<k<<" "; cout<<endl;
    }

    vector<int> r3= Consensus(r1,r2);

    if((count(r3.begin(),r3.end(),0)==0 && r3.size()>0) || (r3.size()==13 && count(r3.begin(),r3.end(),0)==1))
    {
      usable_reads++;
      if(code_readout.find(r3)==code_readout.end()){std::vector<int>I={i}; code_readout[r3]=I;}
      else code_readout[r3].push_back(i);

    }

    cout<<"rating "<<usable_reads/(tot_reads)<<endl;

  }

  std::vector<int> output_code_sizes; output_code_sizes.reserve(code_readout.size());
  std::vector<string> output_code_readout; output_code_readout.reserve(code_readout.size());
  std::vector<int> output_code_counts; output_code_counts.reserve(code_readout.size());
  std::vector<int> output_code_firstread; output_code_firstread.reserve(code_readout.size());
  std::vector<std::vector<int> > output_code_readids;
  if(full) output_code_readids.reserve(code_readout.size());
  for(auto c : code_readout){


    output_code_readout.push_back("");
    for(int i = 0; i < c.first.size(); ++i){ // we suppress start and end
      output_code_readout.back() += std::regex_replace(to_string(c.first[i]),std::regex("^0"), "?");
      // remove trailing ' ' - very important when converting back to integer form
      if(i < c.first.size()-1) output_code_readout.back() += " ";

    }

    cout<<"code "<<output_code_readout.back()<<endl;

    output_code_counts.push_back(c.second.size());
    output_code_firstread.push_back(c.second.front());
    output_code_sizes.push_back(c.first.size());
    if(full) output_code_readids.push_back(c.second);
  }

  Rcpp::DataFrame output_df = Rcpp::DataFrame::create(Named("count") = wrap(output_code_counts),
                                                      Named("firstread") = wrap(output_code_firstread),
                                                      Named("code") = wrap(output_code_readout),
                                                      Named("size") = wrap(output_code_sizes),
                                                      Named("stringsAsFactors") = false);
  Rcpp::S4 decode_output("decode_output");
  decode_output.slot("data") = output_df;
  decode_output.slot("read_ids") = output_code_readids;
  decode_output.slot("saturation") = saturation;
  Rcpp::S4 output("loxcode_sample");

  output.slot("decode") = decode_output;
  output.slot("name") = name;
  output.slot("meta") = meta;
  output.slot("files") = r;

  Rcpp::List decode_stats;
  decode_stats["tot_reads"] = tot_reads;
  decode_stats["too_small"] = too_small_reads;
  //decode_stats["missing_start"] =  reads_missing_start;
  //decode_stats["multi_start"] = reads_multi_start;
  //decode_stats["missing_end"] = reads_missing_end;
  //decode_stats["multi_end"] = reads_multi_end;
  decode_stats["consensus_filtered"] = usable_reads;
  output.slot("decode_stats") = decode_stats;

  output.slot("consensus_filtered_data") = reads_consensus_filtered_data;
  fileR1.close();fileR2.close();

  return output;
}
