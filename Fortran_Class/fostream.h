//fostream.h
#ifndef FOSTREAM_H
#define FOSTREAM_H
// declaration of a SIMPLE iostream-like class
// that is using FORTRAN to implement I/O
// This is useful when doing file I/O in a
// mixed C++/F77 program
#include <stdlib.h>
#include "fortran.h"

/*extern "C" __declspec(dllimport) void _stdcall F77CLS();
extern "C" __declspec(dllimport) void _stdcall F77OUT(CHARACTER STRING);
extern "C" __declspec(dllimport) void _stdcall F77OPN(const INTEGER &IU,CHARACTER NAME);
*/
//------------------------------------------------------------------------------------
SUBROUTINE F77CLS();
SUBROUTINE F77OUT(INTEGER &IUNIT, CHARACTER STRING);
SUBROUTINE F77OPN(const INTEGER &IU,CHARACTER NAME,CHARACTER MODE,int &erro);
SUBROUTINE F77REWIND (const long int &IUNIT,int &erro);
//------------------------------------------------------------------------------------
class fostream {
public:
	int erro,iunit;
    char *mode;
    bool open_mode;
  
    fostream(int funit=6, char* filename=NULL,char *status=NULL);
	 
	int open(int funit, char* filename,char *status);

    int rewind();    

    fostream& operator <<(const char  ch);
    fostream& operator <<(const char* txt);
    ~fostream();
    
};          
//------------------------------------------------------------------------------------
fostream::fostream(int funit, char* filename,char *status)
{
   erro = 1;          //start with 1 to indicate no error
   mode = "UNKNOWN";  //initial mode is unknown
   open_mode = false;     //indicate that the file is not open
   open(funit,filename,mode);
}
//------------------------------------------------------------------------------------
int fostream::rewind()
{
    F77REWIND (iunit,erro);
    return(erro);
}
//------------------------------------------------------------------------------------
int fostream::open(int funit, char* filename,char *status)
{ iunit = funit;
  mode = status;
 if(filename)
  { F77OPN(iunit,CHARACTER(filename),mode,erro);
     if(erro) open_mode = true;
  }
 else
    F77OPN(iunit,CHARACTER("stdout.txt"),mode,erro);



return (erro);

}
//------------------------------------------------------------------------------------
fostream::~fostream()
{
  F77CLS();
}
//------------------------------------------------------------------------------------
fostream& fostream::operator <<(const char ch)
{
   char str[2] = { ch,'\0'};
   F77OUT(iunit, CHARACTER(str));
   return *this;
}
//------------------------------------------------------------------------------------
fostream& fostream::operator <<(const char* txt)
{
   F77OUT(iunit, CHARACTER((char*)txt));
   return *this;
}


#endif

