
//Compute approximate PAM distance
//Take 4 values."_arrid1/_arrid2" AAid sequence."_distype" distance type:'S'simple.'P'poisson.'G'gamma."_gma" gamma value
//Output 4 values."dist"distance(PAM/100)."var" variance."smlt" similarity."BKMATX" Needleman global alignment matrix
void AprxDist(sarray<int> & _arrid1,sarray<int> & _arrid2,double & dist,double & var,double & smlt,const char & _distype,darray<int> & BKMATX,const int &_gma=1)
{
	int match,len;
	NWalign(_arrid1,_arrid2,BKMATX);
    match=BKMATX(0,0);
    len=_arrid1.size()+_arrid2.size();
	if(_distype=='S' || _distype=='s')
	{
		dist=1.0-(double)match*2/(double)len;
		smlt=1.0-dist;
		if(dist==0)
			dist=0.01;
        var=dist*(1.0-dist)/(double)len;
	}
	else if(_distype=='P' || _distype=='p')
	{
		double p;
		p=1.0-(double)match*2/(double)len;
		smlt=1.0-p;
		dist=-1.0*log(1.0-p);
		if(dist==0)
			dist=0.01;
		var=p/((1.0-p)*(double)len);
	}
	else if(_distype=='G' || _distype=='g')
	{
		double p;
		p=1.0-(double)match*2/(double)len;
		smlt=1.0-p;
        dist=_gma*(pow(1.0-p,-1.0/_gma)-1.0);
		if(dist==0)
			dist=0.01;
		var=p*pow(1.0-p,-1.0-2.0/_gma)/(double)len;
	}
    return;
}

HANDLE READY;                               //global var for thread use
//This structure pass parameter from AprxDistMatx to thread function THaprxdist
struct getaprxdiststr
{
	vector<sarray<int>> * _multiarrid;              //pointer to sequence id vec
	darray<double> * _distmatx;                     //pointer to distance matrix
	darray<double> * _varmatx;                      //pointer to variance matrix
	darray<int> * _thrdtaskmatx;                    //pointer to thread task matrix,this matrix makes the task assignment easy
	char _distype;                                  //distance estimation type
	int _gma;                                       //gamma
    int _thrdid;                                    //thread id
};

//Menu:[AprxDistMatx]-[THaprxdist]
//This child thread allocates aprox distance calculation task to all CPUs
//Take pointer to getaprxdiststr structure
//Update "_distmatx" and "_varmatx" in the imported structure
unsigned __stdcall THaprxdist(void * param)
{
	getaprxdiststr *datastr=static_cast<getaprxdiststr *>(param);
	int thrdid=datastr->_thrdid;                                           //assign "_thrdid" to local var
    SetEvent(READY);                                                       //inform the main thread to proceed
	int start,end,size,quantum,I,J;                                        //"start""end"thread task start and end position."quantum" aliquot of task
	double dist,var,smlt;
	size=datastr->_thrdtaskmatx->getcnum();
	quantum=(int)((double)size/CPUNUM);
	if(quantum<1)                                                          //if CPU number is larger than task number,set "quantum" to 1,otherwise it is 0.
		quantum=1;
	start=thrdid*quantum;
	if(thrdid!=CPUNUM-1)
        end=(thrdid+1)*quantum-1;
	else
		end=size-1;
	darray<int> *thrdtaskmatx=datastr->_thrdtaskmatx;
	vector<sarray<int>> *multiarrid=datastr->_multiarrid;
	darray<double> *distmatx=datastr->_distmatx;
	darray<double> *varmatx=datastr->_varmatx;
	darray<int> BKMATX;
    for(int iter=start;iter<=end;++iter)
	{
		if(thrdid==CPUNUM-1)                                               //control output format
			PctMarker(iter-start,end-start,2);
        I=thrdtaskmatx->operator ()(0,iter);                               //the id of seq A
		J=thrdtaskmatx->operator ()(1,iter);                               //the id of seq B
	    AprxDist(multiarrid->operator [](I),multiarrid->operator [](J),dist,var,smlt,datastr->_distype,BKMATX,datastr->_gma);
	    distmatx->operator ()(I,J)=distmatx->operator ()(J,I)=dist;
	    varmatx->operator ()(I,J)=varmatx->operator ()(J,I)=var;
	}
	if(thrdid==CPUNUM-1)                                                   //control output format
	   cout<<endl;
	return 0;
}
//Compute approximate PAM distance for multiple sequences
//Take 3 values."_multiarrid" AAid sequences."_distype" distance type (see AprxDist)."_gma" gamma value
//Output 2 values."distmatx" symmetrical distance matrix."varmatx" symmetrical variance matrix
void AprxDistMatx(vector<sarray<int>> & _multiarrid,darray<double> & distmatx,darray<double> & varmatx,const char & _distype,const int & _gma=1)
{
	int rownum=_multiarrid.size(),iter=0;
	distmatx.fast_resize(rownum,rownum);
	varmatx.fast_resize(rownum,rownum);
	darray<int> thrdtaskmatx(2,(int)((double)rownum*(double)(rownum-1)/2.0));       //this var holds all combinations that needs to be done
	for(int i=0;i<rownum;++i)                                                       //it is specifically designed for multi thread function
	{
		distmatx(i,i)=0.0;
		varmatx(i,i)=0.0;
		for(int j=i+1;j<rownum;++j)
		{
			thrdtaskmatx(0,iter)=i;
            thrdtaskmatx(1,iter)=j;
			iter+=1;
		}
	}
	getaprxdiststr datastr;                                                          //transfer necessary var to structure before enter child thread function
	datastr._multiarrid=&_multiarrid;
	datastr._distmatx=&distmatx;
	datastr._varmatx=&varmatx;
	datastr._thrdtaskmatx=&thrdtaskmatx;
	datastr._distype=_distype;
	datastr._gma=_gma;
	READY=CreateEvent(NULL,FALSE,FALSE,NULL);                                        //global event,for thread synchronization 
	sarray<HANDLE> handlearr;                                                        //handle array holds handles of all child threads
	for(int i=0;i<CPUNUM;++i)                                                        //create "CPUNUM" threads
	{
		datastr._thrdid=i;                                                           //where child threads get their id
		handlearr.pushback((HANDLE)_beginthreadex(NULL,0,&THaprxdist,&datastr,NULL,NULL));    //create child thread and store handle in "handlearr"
		WaitForSingleObject(READY,INFINITE);                                         //wait until child thread signaled to proceed safely
	}
	WaitForMultipleObjects(CPUNUM,handlearr.getarr(),TRUE,INFINITE);                 //wait until all threads doen their work
}
//This structure pass parameter from BysDistMatx to thread function THbysdist
struct getbysdiststr
{
	vector<darray<double>> *_scoredb;             //pointer to score database      
	vector<sarray<int>> *_multiarrid;             //pointer to seq id matrix
	darray<double> *_distmatx;                    //pointer to distance matrix
	darray<double> *_varmatx;                     //pointer to variance matrix
	darray<int> *_thrdtaskmatx;                   //pointer to thread task matrix,this matrix makes the task assignment easy
	int _wlen;                                    //word length
	int _gappnt;                                  //gap penalty
	int _step;                                    //bayesian step
	int _acculvl;                                 //bayesian accuracy level
	int _thrdid;                                  //thread id
	char _aligntype;                              //alignment type
};
//Menu:[BysDistMatx]-[THbysdist]
//This child thread allocates bayesian distance calculation task to all CPUs
//Take pointer to getbysdiststr structure
//Update "_distmatx" and "_varmatx" in the imported structure
unsigned __stdcall THbysdist(void *param)
{
	getbysdiststr *datastr=static_cast<getbysdiststr*>(param);               //notations are same as THaproxdist child thread
	int thrdid=datastr->_thrdid;
	SetEvent(READY);
    int start,end,size,quantum,I,J;
	double dist,var;
	size=datastr->_thrdtaskmatx->getcnum();
	quantum=(int)((double)size/CPUNUM);
	if(quantum<1)
		quantum=1;
	start=thrdid*quantum;
	if(thrdid!=CPUNUM-1)
        end=(thrdid+1)*quantum-1;
	else
		end=size-1;
	vector<darray<double>> *scoredb=datastr->_scoredb;
	vector<sarray<int>> *multiarrid=datastr->_multiarrid;
	darray<double> *distmatx=datastr->_distmatx;
	darray<double> *varmatx=datastr->_varmatx;
	darray<int> *thrdtaskmatx=datastr->_thrdtaskmatx;
    for(int iter=start;iter<=end;++iter)
	{
		if(thrdid==CPUNUM-1)
			PctMarker(iter-start,end-start,2);
        I=thrdtaskmatx->operator ()(0,iter);
		J=thrdtaskmatx->operator ()(1,iter);
		BayesianDist(datastr->_aligntype,*scoredb,multiarrid->operator [](I),multiarrid->operator [](J),dist,var,datastr->_wlen,datastr->_gappnt,datastr->_step,datastr->_acculvl);
		distmatx->operator ()(I,J)=distmatx->operator ()(J,I)=dist;
		varmatx->operator ()(I,J)=varmatx->operator ()(J,I)=var;
	}
	if(thrdid==CPUNUM-1)
		cout<<endl;
	return 0;
}
//Compute bayesian exact PAM distance for multiple sequences
//Take 7 values."_multiarrid" AAid sequences."_scoredb"scoring database(usually mean Dayhoff PAM database)."_aligntype" 'N'for needleman global alignment;'B'for BLAST local alignment."_wlen" word length(for BLAST)."_step"/"_acculvl"(see BayesianDist)
//Output 2 values."distmatx""varmatx"
void BysDistMatx(vector<sarray<int>> & _multiarrid,darray<double> & distmatx,darray<double> & varmatx,vector<darray<double>> & _scoredb,const char & _aligntype,const int & _wlen=3,const int & _gappnt=-3,const int & _step=3,const int & _acculvl=4)
{
	clock_t start,end;                                                       //notations are same as AprxDistMatx
	start=clock();
	BYSINPROGRESS=true;
	int rownum=_multiarrid.size(),iter=0;
	distmatx.fast_resize(rownum,rownum);
	varmatx.fast_resize(rownum,rownum);
	darray<int> thrdtaskmatx(2,(int)((double)rownum*(double)(rownum-1)/2.0));   
	for(int i=0;i<rownum;++i)
	{
		distmatx(i,i)=0.0;
		varmatx(i,i)=0.0;
		for(int j=i+1;j<rownum;++j)
		{
			thrdtaskmatx(0,iter)=i;
            thrdtaskmatx(1,iter)=j;
			iter+=1;
		}
	}
	getbysdiststr datastr;
	datastr._scoredb=&_scoredb;
	datastr._multiarrid=&_multiarrid;
	datastr._distmatx=&distmatx;
	datastr._varmatx=&varmatx;
	datastr._thrdtaskmatx=&thrdtaskmatx;
	datastr._acculvl=_acculvl;
	datastr._aligntype=_aligntype;
	datastr._gappnt=_gappnt;
	datastr._step=_step;
	datastr._wlen=_wlen;
	READY=CreateEvent(NULL,FALSE,FALSE,NULL);
	sarray<HANDLE> handlearr;
	for(int i=0;i<CPUNUM;++i)
	{
		datastr._thrdid=i;
		handlearr.pushback((HANDLE)_beginthreadex(NULL,0,&THbysdist,&datastr,NULL,NULL));
		WaitForSingleObject(READY,INFINITE);
	}
	WaitForMultipleObjects(CPUNUM,handlearr.getarr(),TRUE,INFINITE);
	BYSINPROGRESS=false;
	end=clock();
    /////////////////////////////////////////////////////////////////write log
	char time[10],date[10];
	_strtime(time);_strdate(date);
	ofstream outlog("c:\\ProFunSite.log",ios_base::app);
    outlog<<date<<" "<<time<<" Process time "<<end-start<<" minisec"<<endl;
	outlog<<"Create bayesian exact pairwise distance matrix for "<<rownum<<" sequences  ";
	if(_aligntype=='N' || _aligntype=='n')
		outlog<<"Align type:Needle global alignment  Gap penalty:"<<_gappnt;
	else
		outlog<<"Align type:BLAST local alignment  Word length:"<<_wlen;
	outlog<<"  Step:"<<_step<<"  Accuracy level:"<<_acculvl<<endl<<endl;
}



