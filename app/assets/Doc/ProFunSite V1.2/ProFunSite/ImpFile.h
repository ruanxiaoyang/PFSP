

struct inf                                                                  //Structure for all format.
{
	sarray<char> _gi;
	sarray<char> _ginum;
    sarray<char> _format;
	sarray<char> _accession;
	sarray<char> _locus;
};
//////////////////////////////////////////////////////////////////////////////standard fasta format
struct proformat                                                          //structure hold all information of sequence file
{
	vector<sarray<int>> _multiseqid;
	vector<inf> _infvec;
	int _seqnum;
	char _format;
};
//input fasta format protein sequence
//take 1 value."_path" the directory of profile
//update 2 value."_infvec" information vector of seqs."_multiseq" vector store seqs
int Fasta(char * _path,proformat & seqsandinf)
{
	clock_t start,end;
	start=clock();
	seqsandinf._format='F';
	inf proinf;
	char aa;
	seqsandinf._seqnum=0;
	sarray<int> aaseqid;
	ifstream input(_path);
	if(!input)
		return 0;                                                             //return 0 if the profile was not found
    while(input.get(aa) && aa=='>')
	{
        while(input.get(aa))
		{
			if(aa!='|')
			proinf._gi.pushback(aa);
			else
				break;
		}
		while(input.get(aa))
		{
			if(aa!='|')
			proinf._ginum.pushback(aa);
			else
				break;
		}
        while(input.get(aa))
		{
			if(aa!='|')
			proinf._format.pushback(aa);
			else
				break;
		}
		while(input.get(aa))
		{
			if(aa!='|')
			proinf._accession.pushback(aa);
			else
				break;
		}
		while(input.get(aa))
		{
			if(aa!='\n')
			proinf._locus.pushback(aa);
			else
				break;
		}
		seqsandinf._infvec.push_back(proinf);
		proinf._gi.clear();proinf._ginum.clear();proinf._format.clear();proinf._accession.clear();proinf._locus.clear();
		while(input.get(aa))
		{
			if(aa=='>')
				break;
			if((int)aa>=65 && (int)aa<=90 && aa!='X' && aa!='B' && aa!='Z')
				aaseqid.pushback(aaid(aa));
		}
		seqsandinf._multiseqid.push_back(aaseqid);
		seqsandinf._seqnum+=1;
		aaseqid.clear();
		input.unget();
	}
	end=clock();
	char time[10],date[10];
	_strtime(time);_strdate(date);
	ofstream outlog("c:\\ProFunSite.log",ios_base::app);
    outlog<<date<<" "<<time<<" Process time "<<end-start<<" minisec"<<endl;
	outlog<<seqsandinf._multiseqid.size()<<" sequences uccessfully imported"<<endl<<endl; 
	input.close();
	return 1;
}

/////////////////////////////////////////////////////////////////////////////custom format,same structure but different name

int Custom(char *_path,proformat & seqsandinf)
{
	clock_t start,end;
	start=clock();
	seqsandinf._format='C';
	char aa;
	seqsandinf._seqnum=0;
	sarray<int> aaseqid;
	inf proinf;
	ifstream input(_path);
	if(!input)
		return 0;                                                             //return 0 if the profile was not found
    while(input.get(aa) && aa=='>')
	{
		while(input.get(aa))
		{
			if(aa!='\n')
				proinf._gi.pushback(aa);
			else
				break;
		}
		seqsandinf._infvec.push_back(proinf);
		proinf._gi.clear();
		while(input.get(aa))
		{
			if(aa=='>')
				break;
			if((int)aa>=65 && (int)aa<=90 && aa!='X' && aa!='B' && aa!='Z')
				aaseqid.pushback(aaid(aa));
		}
		seqsandinf._multiseqid.push_back(aaseqid);
		seqsandinf._seqnum+=1;
		aaseqid.clear();
		input.unget();
	}
	end=clock();
	char time[10],date[10];
	_strtime(time);_strdate(date);
	ofstream outlog("c:\\ProFunSite.log",ios_base::app);
    outlog<<date<<" "<<time<<" Process time "<<end-start<<" minisec"<<endl;
	outlog<<seqsandinf._multiseqid.size()<<" sequences uccessfully imported"<<endl<<endl; 
	input.close();
	return 1;
}

