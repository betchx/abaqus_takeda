#include <iostream>
#include <fstream>
#include <windows.h>
#include <boost/algorithm/string.hpp>

using namespace std;

#ifndef MAX_PATH
#define MAX_PATH 1024
#endif

int main()
{
  DWORD len;
  len = GetCurrentDirectory(0,NULL);
  if(len == 0) {
    cerr << "Could not get Current Directory" << endl;
    return -1;
  }
  TCHAR* dirname = new TCHAR[len];
  if(!dirname){
    cerr << "Alloc Error" << endl;
    return -2;
  }
  if(GetCurrentDirectory(len,dirname)==0){
    delete[] dirname;
    cerr << "Could not get Current Directory" << endl;
    return -3;
  }
  // Substitute '\\' to '/'
  for(TCHAR* ch = dirname; *ch != '\0'; ++ch){
    if(*ch == '\\') *ch = '/';
  }
  ofstream ofs("abaqus_v6.env");
  ofs << "ask_delete=OFF" << endl;
  ofs << endl;
  ofs << "link_sl='cmd /c \""
    << "LINK /nologo"
    << " /MANIFEST"
    << " /INCREMENTAL:NO"
    << " /subsystem:console"
    << " /machine:X86"
    << " /LIBPATH:" << dirname;
  delete[] dirname;
  ofs << " /DEFAULTLIB:uels.lib"
      << " /NODEFAULTLIB:LIBC.LIB"
      << " /NODEFAULTLIB:LIBCMT.LIB"
      << " /NODEFAULTLIB:LIBIFCOREMT.LIB"
      << " /NODEFAULTLIB:libmmt.lib"
      << " /DEFAULTLIB:OLDNAMES.LIB"
      << " /DEFAULTLIB:LIBIFCOREMD.LIB"
      << " /DEFAULTLIB:LIBIFPORTMD.LIB"
      << " /DEFAULTLIB:LIBMMD.LIB"
      << " /DEFAULTLIB:MSVCRT.LIB"
      << " /DEFAULTLIB:kernel32.lib"
      << " /DEFAULTLIB:user32.lib"
      << " /DEFAULTLIB:advapi32.lib"
  //    << " /LIBPATH:C:\\Intel\\Compiler\\11.0\066\\fortran\\mkl\\ia32\\lib"
  //    << " /DEFAULTLIB:mkl_intel_s.lib"
  //    << " /DEFAULTLIB:mkl_intel_thread.lib"
  //    << " /DEFAULTLIB:mkl_sequential.lib"
  //    << " /DEFAULTLIB:mkl_core.lib"
      << " /DEFAULTLIB:libiomp5md.lib"
      << " /FIXED:NO"
      << " /dll"
      << " /def:%E"
      << " /out:%U %F %A %B"
      << " && mt /manifest %U.manifest /outputresource:%U;2"
      << " && del %U.manifest\"'" << endl;
  ofs.close();
  return 0;
}


