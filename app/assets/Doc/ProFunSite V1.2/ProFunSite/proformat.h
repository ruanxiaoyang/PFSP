


struct inf                                     //structure specific for fasta format
{
	sarray<char> _gi;
	sarray<char> _ginum;
    sarray<char> _format;
	sarray<char> _accession;
	sarray<char> _locus;
};

//input fasta format protein sequence
//take 1 value."_path" the directory of profile
//update 2 value."_infvec" information vector of seqs."_multiarray" vector store seqs
int Fasta(char * _path,vector<inf> & _infvec,vector<chararray> & _multiarray)
{
	inf proinf;
	char aa;
	sarray<char> aaseq;
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
		_infvec.push_back(proinf);
		proinf._gi.clear();proinf._ginum.clear();proinf._format.clear();proinf._accession.clear();proinf._locus.clear();
		while(input.get(aa))
		{
			if(aa=='>')
				break;
			if(aa!='\n' && aa!=' ' && aa!='-')
				aaseq.pushback(aa);
		}
		_multiarray.push_back(aaseq);
		aaseq.clear();
		input.unget();
	}
	return 1;
}