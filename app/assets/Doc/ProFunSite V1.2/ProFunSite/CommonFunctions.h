
//Calculate log standard binomial probability
//Take 3 values."_N" total number."_n"number of occurance."_p"probability
//Return 1 value.The log probability of "_n"occurances out of "_N" times when total probability is "_p".To obtain p,use exp(LogBinoP(...))
inline double LogBinoP(const int & _N,const int & _n,const double & _p)
{
	if(_p==0)
		return 0.0;
	double nt=0,Nt=0;
	for(int i=_N;i>_n;--i)
		Nt+=log((double)i);
	for(int i=_N-_n;i>=1;--i)
		nt+=log((double)i);
	return  (Nt-nt)+_n*log(_p)+(_N-_n)*log((double)(1-_p));
}

//Control format.Output number mark for sequence.Spaced by 10 AA/nucleotide.
//Take 2 values."_size"length of sequence."_os"output file stream
inline void LandMarker(const int & _size,ostream & _os)
{
	int hang=0;
	for(int i=0;i<_size;++i)
	{
		if(i>=10 && i%10==0)
		{
			_os<<i;
		    hang=(int)(log((double)i)/log(10.0));
			continue;
		}
		if(hang==0)
            _os<<" ";
		else if(hang>0)
			hang-=1;
	}
}
//Create multi-level directory
//Take 1 value."_path" character array (such as "C:\\TEMP")
void MKDIR(char* _path)
{
	bool enter=false;
	string dir,path(_path);                           //transfer char array to string
	for(int i=0;i<path.size();++i)
	{
		if(path[i]!='\\' && enter==false)
			dir.push_back(path[i]);
		if(path[i]=='\\' && enter==false)
		{
			dir.push_back(path[i]);
			enter=true;
			continue;
		}
		if(enter==true)
		{
			if(path[i]=='\\')
			{
				_mkdir(dir.c_str());                  //_mkdir is a C function,which takes char*(directory).it can only create single-level directory
				dir.push_back(path[i]);
			}
			else if(i==path.size()-1)
			{
				dir.push_back(path[i]);
		        _mkdir(dir.c_str());
			}
			else
				dir.push_back(path[i]);
		}
	}
}

//Control format.Output "+" mark for sequences.Spaced by 10 AA/nucleotide
//Take 2 values."_size"length of sequence."_os"output file stream
inline void MileStone(const int & _size,ostream & _os)
{
	for(int i=0;i<_size;++i)
	{
		if(i%10==0)
			_os<<"+";
		else
			_os<<"-";
	}
}

//Control format.Same as progressbar.
//Take 2 values."_i"usually the start number of iteration."_size"usually the end number of iteration."_acry" round to _acry digits after decimal point
inline void PctMarker(const int & _i,const int & _size,const int & _acry=0)
{
	if(_i>_size)
		return;
	double per;
	if(_size==0)
		per=1;
	else 
		per=(double)_i/_size;
	int tmp=(int)(per*pow(10.0,(double)(2+_acry)))*10;
	cout.flags(ios::fixed);
	cout.width(3+_acry);
	cout.precision(_acry);
	if(tmp%10==0)
	{
			cout<<tmp/pow(10.0,(double)(_acry+1))<<"%";
			if(_i!=_size)
			{
				for(int i=0;i<4+_acry;++i)
					cout<<"\b";
			}
	}
}
//Control file name,change int type to string type
inline string inttostr(const int & temp)
{
	stringstream stream;
	string g;
	stream<<temp;
	stream>>g;
	return g;
}

//compute length of sequence after alignment(gap+original length)
//take 3 values."_matrix" standard 0 or 1 matrix."_rownum"/"_colnum" row and column of "_matrix"
//return 1 value.length of sequence after alignment
int Length(darray<int> & _matrix,const int & _rownum,const int & _colnum)
{
	int previ=-1,prevj=-1,gap=0;
	for(int i=0;i<_rownum;++i)
	{
		for(int j=prevj+1;j<_colnum;++j)
		{
			if(_matrix(i,j)==1)
			{
				gap+=i-previ-1;                                //compute row gap 
				previ=i;prevj=j;
				if(prevj==_colnum-1)
				{
					gap+=_rownum-1-i;
					return gap+_colnum;
				}
				break;
			}
		}
	}
	return gap+_colnum;
}

//Break string obtain from getline() function into blank ' 'seperated words and store these words in vector<sarray<char>>
//Take 1 value. "_tmpstr" string type
//Return 1 value.vector<sarray<char>> containing word of "_tempstr"
vector<sarray<char>> _Breakstring(string & _tmpstr)
{
	vector<sarray<char>> wordvec;
	sarray<char> word;
	for(int i=0;i<_tmpstr.size();++i)
	{
		if(_tmpstr[i]!=' ')
		{
			word.pushback(_tmpstr[i]);
			continue;
		}
		if(_tmpstr[i]==' ' && word.size()==0)
			continue;
		if(_tmpstr[i]==' '&& word.size()>0)
		{
			wordvec.push_back(word);
			word.clear();
		}
	}
	if(word.size()>0)
		wordvec.push_back(word);
	return wordvec;
}
//Copy file from srcFile to destFile
//"srcFile/destFile"are char array.Example Copy("blocks.dat","c:\\windows\\blocks.dat");
void Copy(const char* srcFile, const char* destFile)
{
    ifstream in(srcFile);
    ofstream out(destFile);
    out<<in.rdbuf();
    out.close();
    in.close();
}
//Convert AA id seqs to AA char seqs
//Take 1 value."_MatxInt" AA id seqs
//Output 1 value."MatxChar" AA char seqs
void MatxItoC(darray<int> & _MatxInt,darray<char> & MatxChar)
{
	MatxChar.fast_resize(_MatxInt.getrnum(),_MatxInt.getcnum());
	for(int i=0;i<_MatxInt.getrnum();++i)
	{
		for(int j=0;j<_MatxInt.getcnum();++j)
			MatxChar(i,j)=(_MatxInt(i,j)>=0 ? idaa(_MatxInt(i,j)):'-');
	}
}
//Convert AA char seqs to AA id seqs
void MatxCtoI(darray<char> & _MatxChar,darray<int> & MatxInt)
{
	MatxInt.fast_resize(_MatxChar.getrnum(),_MatxChar.getcnum());
	for(int i=0;i<_MatxChar.getrnum();++i)
	{
		for(int j=0;j<_MatxChar.getcnum();++j)
			MatxInt(i,j)=aaid(_MatxChar(i,j));
	}
}
//This structure holds file read from _ProFunSit_t_s_thld
struct _TSSTR
{
	char _DorB;                         //Dayhoff or BLOSUM scoring matrix
	int _scmatxid;                      //id of scoring matrix (PAM250)(BLOSUM60)
	int _wlen;                          //word length
	int _seqlen;                        //sequence length
	int _cutoffper;                     //cutoff percentile,usually 95
	double _tthld;
	double _sthld;
};

vector<_TSSTR> TSSTRVEC;                //global var holds file read from _ProFunSit_t_s_thld
//This function read _ProFunSit_t_s_thld file into memory
//Take 1 value."_imppath" the file path of _ProFunSit_t_s_thld
//_ProFunSit_t_s_thld file is stored in global var "TSSTRVEC"
void readTSintomemory(string & _imppath)
{
	ifstream read(_imppath.c_str());             //look up existing file for t and s score threshold
	char a;
	int b,c,d,e;
	double tthld,sthld;
	_TSSTR temp;
	if(read)
	{
		while(read>>a)
		{
			read>>b>>c>>d>>e>>tthld>>sthld;
			temp._DorB=a;
			temp._scmatxid=b;
			temp._wlen=c;
			temp._seqlen=d;
			temp._cutoffper=e;
			temp._tthld=tthld;
			temp._sthld=sthld;
			TSSTRVEC.push_back(temp);
		}
	}
}
//This function write "TSSTRVEC" to file _ProFunSit_t_s_thld
void writeTStofile()
{
	ofstream write("c:\\windows\\_ProFunSit_t_s_thld");
	for(int i=0;i<TSSTRVEC.size();++i)
		write<<TSSTRVEC[i]._DorB<<" "<<TSSTRVEC[i]._scmatxid<<" "<<TSSTRVEC[i]._wlen<<" "<<TSSTRVEC[i]._seqlen<<" "<<TSSTRVEC[i]._cutoffper<<" "<<TSSTRVEC[i]._tthld<<" "<<TSSTRVEC[i]._sthld<<endl;
}

//Get CPU number
int core_count()
{
  int count = 1;
  SYSTEM_INFO si;
  GetSystemInfo(&si);
  count = si.dwNumberOfProcessors;
  return count;
}

DWORD WINAPI dspfile(void *path)
{
	string *pathtmp=static_cast<string*>(path);
	system(pathtmp->c_str());
	return NULL;
}

void DSPFILE(string & path)
{
	CreateThread(NULL,0,&dspfile,&path,NULL,NULL);
}

