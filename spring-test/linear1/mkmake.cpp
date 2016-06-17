#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char* argv[])
{
  vector<string> lst;

  string line;
  string target;

  while(cin >> line)
  {
    lst.push_back(line);
  }

  if(lst.size() == 0){
    cerr << "解析対象候補となるinpファイルがみつかりません." << endl;
    return -1;
  }

  if(lst.size() > 1){
    int i = 1;
    int res = 0;
    while(true){
      cerr << "解析対象のInpファイルはどれですか？"<<endl;
      cerr << "  0以下: キャンセル" << endl;
      for(vector<string>::iterator it = lst.begin(); it != lst.end(); it++,i++){
        cerr << "  " << i << ":" << *it << endl;
      }
      cin >> res;
      if(res < 1) return -1;
      if(res <= lst.size()) break;
      cerr << "指定が無効です．" << endl;
      return -1;
    }
    target = lst[res - 1];
  }else{
    target = lst[0];
  }

//  cout << "all: " <<target << " wheels.inp Param1.DAT Param2.DAT Param3.DAT ejrs.for ejrs.ctl FIR_LP_0.01_N105.txt" << endl;
//  cout << "	abaqus user=ejrs interactive job=" << target.substr(0, target.size() - 4) << endl << endl;
  cout << "all: " <<target << " uel.for" << endl;
  cout << "	abq6114 user=uel ";
  if(argc > 0) cout << "interactive ";
  cout << "job=" << target.substr(0, target.size() - 4) << endl << endl;

  return 0;

}



