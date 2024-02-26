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

std::map<string, int> ele_R={{"ACTCCGCA",1},{"TCCAGAATTTGTAT",2},{"ACATCCAC",3},{"AAAGGAATTTCTCC",4},{"ATTTCCTC",5},{"GCCCGAATTTTTTC",6},{"GCTACTGG",7},{"ATGAGAATTTATGG",8},{"AACTAGAA",9},{"TGCAGAATTTCCTC",10},{"CGACACTT",11},{"AACGGAATTTTCAA",12},{"CGTGTTTG",13},{"ACACGAATTCATCC",14},{"GTGCAAATTCGAGA",-15},{"CAAACACG",-13},{"TTGAAAATTCCGTT",-12},{"AAGTGTCG",-11},{"GAGGAAATTCTGCA",-10},{"TTCTAGTT",-9},{"CCATAAATTCTCAT",-8},{"CCAGTAGC",-7},{"GAAAAAATTCGGGC",-6},{"GAGGAAAT",-5},{"GGAGAAATTCCTTT",-4},{"GTGGATGT",-3},{"ATACAAATTCTGGA",-2},{"TGCGGAGT",-1}};

template<typename T>
void pop_front(std::vector<T>& vec)
{
  assert(!vec.empty());
  vec.erase(vec.begin());
}

std::vector<int> Consensus_SP500(vector<int> a, vector<int> b)
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
    //CHANGE:: 12 -> 13
    if(a.size()+b.size()==13 && count(a.begin(),a.end(),0)==0 && count(b.begin(),b.end(),0)==0)
    {
      //CHANGE:: check for zero intersection:
      sort(a.begin(), a.end());
      sort(b.begin(), b.end());
      vector<int> v(a.size() + b.size());
      vector<int>::iterator it, st;
      set_intersection(a.begin(),a.end(),b.begin(),b.end(),v.begin());
      if(v.size()-count(v.begin(),v.end(),0)>0) return empty;
      /////////////////////////////////////////////////////////

      vector<int> full_length(a);
      //full_length.push_back(0);
      full_length.insert(full_length.end(), b.begin(), b.end());
      //cout<<"full length"<<endl; for(auto i : full_length) cout<<i<<" "; cout<<endl;
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

vector<int> align_R(string R){

  EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_LOC, NULL, 0);
  //printf("A");
  vector<string> LoxPs ={R.substr(0,34),R.substr(42,34),R.substr(90,34),R.substr(132,34),R.substr(180,34)};//,R.substr(222,34),R.substr(270,34)};
  //printf("B");
  vector<string> ELE   ={R.substr(34,8),R.substr(76,14),R.substr(124,8),R.substr(166,14),R.substr(214,8)}; //, R.substr(256,14),R.substr(304,8)};
  //printf("C");
  string LoxF="ATAACTTCGTATAATGTATGCTATACGAAGTTAT";
  string LoxR="ATAACTTCGTATAGCATACATTATACGAAGTTAT";

  vector<int> elements;
  for(int k=0; k<LoxPs.size(); k++){
    //printf("D");
    EdlibAlignResult rF,rR;
    rF=edlibAlign(&LoxF[0u],LoxF.size(),&LoxPs[k][0u],LoxPs[k].size(), config);
    //printf("E");
    rR=edlibAlign(&LoxR[0u],LoxR.size(),&LoxPs[k][0u],LoxPs[k].size(), config);
    //printf("F");

    if(std::min(rF.editDistance,rR.editDistance)<10){//likely to be a LoxP site
      //printf("G");
      auto it = ele_R.find(ELE[k]); if(it!=ele_R.end() && it->first.size() == ELE[k].size()) {elements.push_back(it->second); continue;}//exact match
      else{
        bool found=false;
        for(auto e: ele_R) {
          EdlibAlignResult rE=edlibAlign(&ELE[k][0u],ELE[k].size(),&e.first[0u],e.first.size(), config);
          int d=rE.editDistance;
          edlibFreeAlignResult(rE);
          vector<int> exclude={1,3,5,2,8};
          if(d==1 && find(exclude.begin(),exclude.end(),e.second)==exclude.end()) {elements.push_back(e.second); found=true;continue;}//edit distance 1}
        }
        if(found==false) elements.push_back(0);
      }}
    else elements.push_back(0);
    //printf("H");
    edlibFreeAlignResult(rF);edlibFreeAlignResult(rR);
  }
  return elements;
}

void readFASTA_SP500(ifstream &F , std::vector<string> & lines)
{
  getline(F,lines[0]);getline(F,lines[1]);getline(F,lines[2]);getline(F,lines[3]);
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
Rcpp::S4 decode_SP500(std::vector<std::string> files, std::string name, Rcpp::DataFrame meta,
                int min_r1_len, int min_r2_len, bool full, bool sat){
  printf("here");

  double total_reads = 0;
  int reads_consensus_filtered = 0;

  int ok_reads=0;

  map<std::vector<int>,pair<int,int> > count_loxcodes; //loxcode | first read | count

  ifstream fileR1(files[0]); ifstream fileR2(files[1]);
  if(!fileR1.is_open() || !fileR2.is_open()){cout<<"file not found "<<endl;}

  for(int i=0; ;i++){

    //if(i>500000) break;

    // i is the read counter (starts from zero)
    if(fileR1.eof() || fileR2.eof()) break;
    std::vector<string> lines(4);
    readFASTA_SP500(fileR1, lines);
    string R1=lines[1];
    readFASTA_SP500(fileR2, lines);
    string R2=lines[1];
    //parameters specific to SP500
    total_reads++;

    /////////////////trim to start of Lox////////////////////////////////
    int start_R1 = R1.find("TACCGAGCTCGAATTTGCAC"); //20bp before LoxCode R1
    if(start_R1==std::string::npos) continue;
    string R1_trim = R1.substr(start_R1+20,string::npos);

    int start_R2 = R2.find("GCGCCTGGATGAATTCGTGT"); //20bp before LoxCode R2
    if(start_R2==std::string::npos) continue;
    string R2_trim =R2.substr(start_R2+20,string::npos);
    ///////////////////////////////////////////////////

    if(R1_trim.size()<222 || R2_trim.size()<222) continue;
    //printf("********1");
    //printf("%ld\n",R2.length());
    //printf("%ld",R2_trim.length());
    vector<int> r1=  align_R(R1_trim);
    vector<int> r2=  align_R(R2_trim);
    //printf("********2");
    std::reverse(r2.begin(),r2.end()); for(int j=0; j<r2.size(); j++) r2[j]*=-1;
    //printf("********3");
    vector<int> r3= Consensus_SP500(r1,r2);
    //printf("********4");
    set<int> cons_test(r3.begin(),r3.end());
    //no unrecognized elements, size>0 and no repeated elements
    //printf("********5");
    if(count(r3.begin(),r3.end(),0)==0 && r3.size()>0 && cons_test.size()==r3.size()){
      //printf("********6");
      if (count_loxcodes.find(r3) == count_loxcodes.end()) count_loxcodes[r3] = make_pair(i,1);
      else count_loxcodes[r3].second++;
      //printf("********7");
      reads_consensus_filtered++;
    }
    //printf("********8");
  }

////////////////////////////////////////////
  std::vector<int> output_code_sizes; output_code_sizes.reserve(count_loxcodes.size());
  std::vector<string> output_code_elements; output_code_elements.reserve(count_loxcodes.size());
  std::vector<int> output_code_counts; output_code_counts.reserve(count_loxcodes.size());
  std::vector<int> output_code_firstread; output_code_firstread.reserve(count_loxcodes.size());

  for(auto c : count_loxcodes){
    output_code_elements.push_back("");
    for(int i = 0; i < c.first.size(); ++i){ // we suppress start and end
      output_code_elements.back() += std::regex_replace(to_string(c.first[i]),std::regex("^0"), "?");
      // remove trailing ' ' - very important when converting back to integer form
      if(i < c.first.size()-1) output_code_elements.back() += " ";
    }
    output_code_sizes.push_back(c.first.size());
    output_code_firstread.push_back(c.second.first);
    output_code_counts.push_back(c.second.second);
  }


  Rcpp::DataFrame output_df = Rcpp::DataFrame::create(Named("count") = wrap(output_code_counts),
                                                      Named("firstread") = wrap(output_code_firstread),
                                                      Named("code") = wrap(output_code_elements),
                                                      Named("size") = wrap(output_code_sizes),
                                                      Named("stringsAsFactors") = false);
  Rcpp::S4 decode_output("decode_output");
  decode_output.slot("data") = output_df;

  Rcpp::S4 output("loxcode_sample");

  output.slot("decode") = decode_output;
  output.slot("name") = name;
  output.slot("meta") = meta;
  output.slot("files") = files;

  Rcpp::List decode_stats;
  decode_stats["tot_reads"] = total_reads;
  decode_stats["consensus_filtered"] = reads_consensus_filtered;
  output.slot("decode_stats") = decode_stats;

  return output;
}
