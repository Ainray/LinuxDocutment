#include "stdafx.h"
#include "shiwen_rae.h"
#ifdef _DEBUG
#include<time.h>
static int porosity_mul=0;
#endif
// 默认pvt曲线
#ifdef _DOUBLE
static REAL pvtp[]={23.785456,25.738695,29.471909,31.56692,33.583164,
	37.284875,45.995704,56.549523,70.285222,92.337957,109.192549,124.393186,
	141.184773,151.281775,169.711565,189.811062,209.879055,233.317963,250.046535,
	273.45394,290.166771,311.888709,326.947585,353.678651,370.391472,395.468583,
	413.835368,432.202143,453.90834,477.299994,495.139269,562.447418,582.185878,1267.209132,1293.757361};
static REAL pvtv[]={0.055131,0.051737,0.047026,0.041935,0.037787,0.033453,0.029127,
	0.022728,0.018221,0.014103,0.012242,0.01019,0.009084,0.008156,0.007429,0.006705,
	0.006359,0.00564,0.005288,0.004947,0.004783,0.004628,0.004274,0.004126,0.003963,
	0.00368,0.003652,0.003624,0.0035,0.0033,0.00294,0.0028,0.00277,0.00221,0.00221};
static int pvtn=sizeof(pvtv)/sizeof(REAL);
#else
static REAL pvtp[]={23.785456f,25.738695f,29.471909f,31.56692f,33.583164f,
	37.284875f,45.995704f,56.549523f,70.285222f,92.337957f,109.192549f,124.393186f,
	141.184773f,151.281775f,169.711565f,189.811062f,209.879055f,233.317963f,250.046535f,
	273.45394f,290.166771f,311.888709f,326.947585f,353.678651f,370.391472f,395.468583f,
	413.835368f,432.202143f,453.90834f,477.299994f,495.139269f,562.447418f,582.185878f,1267.209132f,1293.757361f};
static REAL pvtv[]={0.055131f,0.051737f,0.047026f,0.041935f,0.037787f,0.033453f,0.029127f,
	0.022728f,0.018221f,0.014103f,0.012242f,0.01019f,0.009084f,0.008156f,0.007429f,0.006705f,
	0.006359f,0.00564f,0.005288f,0.004947f,0.004783f,0.004628f,0.004274f,0.004126f,0.003963f,
	0.00368f,0.003652f,0.003624f,0.0035f,0.0033f,0.00294f,0.0028f,0.00277f,0.00221f,0.00221f};
static int pvtn=sizeof(pvtv)/sizeof(REAL);

#endif

#ifdef _DEBUG
	static REAL maxoilgen=0.0;
#endif
// 通用函数
REAL EXP(REAL x){
#ifdef _DOUBLE
	return exp(x);
#else
	return expf(x);
#endif
}
REAL COS(REAL x){
#ifdef _DOUBLE
	return cos(x);
#else
	return cosf(x);
#endif
}
REAL FABS(REAL x){
#ifdef _DOUBLE
	return fabs(x);
#else
	return fabsf(x);
#endif
}
REAL RANDOM(){
	INT tmp=rand()%100;
	return (REAL) (tmp*1.0/100);
}
REAL RANDOM(REAL a, REAL b){
	return RANDOM()*(b-a)+a;
}
vector<REAL> vectori2r(vector<INT>& x){
	vector<REAL> tmp;
	for(vector<INT>::iterator it=x.begin();it!=x.end();++it)
		tmp.push_back(static_cast<REAL>(*it));
	return tmp;
}
vector<INT> vectorr2i(REAL *x,int n){
	vector<INT> tmp;
	for(int i=0;i<n;++i)
		tmp.push_back(static_cast<INT>(*(x+i)));
	return tmp;
}
CString string2cstring(string str){
	CString tmp(str.c_str());
	return tmp;
}
CString int2cstring(int n){
	CString tmp;
	tmp.Format(_T("%d"),n);
	return tmp;
}

//连通体成员函数
connectorgraph::connectorgraph(){
	_volume_=(REAL)0.0;
	_validvolume_=0.0;
	_oillost_=(REAL)0.0;
	_srcamount_=(REAL) 0.0;
}
connectorgraph::connectorgraph(PCONNECTORGRID pcgd){
	_volume_=pcgd->volume;
	_validvolume_=pcgd->validvolume;
	_oillost_=(REAL) 0.0;
	_srcamount_=(REAL) 0.0;
	_grids_.push_back(pcgd);
}
void connectorgraph::append(PCONNECTORGRID pcgd){
//	vector<connectorgrid>::iterator it=std::upper_bound(_grids_.begin(),_grids_.end(),cgd,DescendingSortByIndex());
	//_grids_.insert(it,cgd);  // insert destroy the iterator, so distance must be calculated before it
	_grids_.push_back(pcgd);
	_volume_+=pcgd->volume;
	_validvolume_+=pcgd->validvolume;
} 
void connectorgraph::createpath(){
	_path_.push_back(vector<INT>());
}
void connectorgraph::appendpath(vector<INT>& path){
	_path_.back().insert(_path_.back().end(),path.begin(),path.end());
}
void connectorgraph::poppath(){
	_path_.pop_back();
}

void connectorgraph::clear(){
	_grids_.swap(vector<PCONNECTORGRID>());
	for(vector<vector<INT>>::iterator it=_path_.begin();it!=_path_.end();++it)
		it->swap(vector<INT>());
	_path_.swap(vector<vector<INT>>());
	delete this;
}

//网格整体类成员函数
shiwen_connector::shiwen_connector(int t){
	_pd_=NULL;
	_type_=t;
}
shiwen_connector::~shiwen_connector(){
	this->clear();
	delete _pd_;
}
int shiwen_connector::initialize(){
	if(_pd_){
		_oilgenpos_=OILGEN_POS;
		_porositypos_=POROSITY_POS;
		_zcoorpos_=ZCOOR_POS;		
		_nx_=this->_pd_->ShiWenFileHead.nXNodeNum-1;
		_ny_=this->_pd_->ShiWenFileHead.nYNodeNum-1;
		//---------------20190303----------------------------------
		rl2cstring();
		_nn_=(this->_pd_->ShiWenFileHead.nXNodeNum-1)*(this->_pd_->ShiWenFileHead.nYNodeNum-1)*(_nz_);
		//---------------------------------------------------------
		_nnx_=this->_pd_->ShiWenFileHead.nXNodeNum;
		_nny_=this->_pd_->ShiWenFileHead.nYNodeNum;
		_nnz_=this->_pd_->ShiWenFileHead.nZNodeNum;
		_nnn_=this->_pd_->ShiWenFileHead.nXNodeNum*this->_pd_->ShiWenFileHead.nYNodeNum*this->_pd_->ShiWenFileHead.nZNodeNum;	
		_toplayer_=0;
		_srclayer_=_nz_-1;
		_ncp_=-1;
		_bottomlayer_=_srclayer_-1;
		_pressurecoe_=SHIWEN_PRESSURE_COE;
		_porositylimit_=SHIWEN_POROSITY;
		_constrainedwater_=SHIWEN_CONTRAINEDWATER;
		_waterden_=WATER_DENSITY;
		_gravityacc_=GRAVITY_ACC;
		_tension_=SHIWEN_TENSION;
		_anglecontact_=ANGLE_CONTACT;	
		 _pathnumberall_=0;
		_surfaltitude_=vector<REAL>(_nx_*_ny_,SHIWEN_SURFACE_ALTITUDE);
		_sandstone_=vector<REAL>(_nnn_,SHIWEN_SANDSTONE_PC);
		_ssstatus_=0;
		_zcoor_.reserve(_nnn_);
		_porosity_.reserve(_nnn_);
		_maxporosity_=0.0;
		_minimaxporosity_=0.0;
		_porosityformat_=1; // % default
		_oilgen_.reserve(_nnn_);
		_oilabundance_=vector<REAL>(_nn_,0);
		_cgds_.reserve(_nn_);
		_gridconnectormap_.cn=vector<INT>(_nn_,-1);
		_gridconnectormap_.pos=vector<INT>(_nn_,-1);
		reinit();
		return TRUE;
	}
	else
		return FALSE;
}
void shiwen_connector::reinit(){
	if(_type_==TYPE_OIL){
		_oilden_=SHIWEN_OIL_DENSITY;
		_oildischargecoe_=SHIWEN_OIL_DISCHARGE_COE;
		_oilgenmin_=SHIWEN_OIL_GEN_INTENSITY_MIN;
		_residualoil_=SHIWEN_RESIDUALOIL;
	}
	else{
		_oilden_=SHIWEN_GASDEN;
		_oildischargecoe_=SHIWEN_GAS_DISCHARGE_COE;
		_oilgenmin_=SHIWEN_GAS_GEN_INTENSITY_MIN;
		_residualoil_=SHIWEN_DISSOLVEDGAS;
		_stoneden_=SHIWEN_STONEDEN;
		_adsorbedgas_=SHIWEN_ADSORBEDGAS;
		_pressureallowederr_=SHIWEN_PERSSURE_ERROR;
		_pvtp_.clear();
		_pvtp_.insert(_pvtp_.begin(),pvtp,pvtp+pvtn);
		_pvtv_.clear();
		_pvtv_.insert(_pvtv_.begin(),pvtv,pvtv+pvtn);
	}
}
void shiwen_connector::flushconnector(){
	for(vector<PCONNECTORGRAPH>::iterator it=_connector_.begin();it!=_connector_.end();++it)
		(*it)->clear();
	vector<PCONNECTORGRAPH>().swap(_connector_);
	// flush buffer
	_cgds_.clear();
	for(int i=0;i<_nn_;++i){
		_gridconnectormap_.cn[i]=-1;
		_gridconnectormap_.pos[i]=-1;
	}
	_minimaxporosity_=0.0;
}
void shiwen_connector::clear(){
	for(vector<PCONNECTORGRAPH>::iterator it=_connector_.begin();it!=_connector_.end();++it)
		(*it)->clear();
	vector<PCONNECTORGRAPH>().swap(_connector_);
	
	vector<CString>().swap(_rlname_);
	vector<INT>().swap(_gridconnectormap_.cn);
	vector<INT>().swap(_gridconnectormap_.pos);
	vector<CONNECTORGRID>().swap(_cgds_);
	vector<REAL>().swap(_oilabundance_);		
	vector<REAL>().swap(_surfaltitude_);
	vector<REAL>().swap(_sandstone_);
	vector<REAL>().swap(_pvtp_);
	vector<REAL>().swap(_pvtv_);
}
int shiwen_connector::read_sproperty(const string &fname,vector<REAL>&pp){
		FILE *pfile;
		REAL tmp;
		int err=fopen_s(&pfile,fname.c_str(),"r");
		if(err) 
			return 0; // error happened
		pp.clear();
		/*	for(int i=0;i<n;++i)
		{
			fscanf_s(pfile,"%f",&tmp);
			pp.push_back(tmp);
		}*/
		int nn=0;
		int c;
		while((c=fscanf_s(pfile,"%f",&tmp))!=EOF){
			pp.push_back(tmp);
			nn++;
		}
		if(c==EOF && nn!=_nnn_){//number not correct
			fclose(pfile);
			return nn;
		}	
		return -1;
}
int shiwen_connector::read_property(const string &fname){
	if(!_pd_){
		_pd_=new iPeraModData();
		if(OnReadShiWenData(fname, *_pd_)){
			initialize();
#ifdef _DEBUG 
			porosity_mul=0;
#endif
			return TRUE;
		}
		else
			return FALSE;
	}
	else
		return FALSE;
}
int shiwen_connector::read_pvt(const string &fname){
	int dcount = 0;
	REAL pvalue;
	REAL vvalue;
	FILE *file;
	_pvtp_.clear();
	_pvtv_.clear();
	fopen_s(&file,fname.c_str(),"r");
	if (file)
	{
		char buffer[256];
		fgets(buffer,256,file);
		while (!feof(file))
		{
			if (fgets(buffer,256,file))
			{
#ifdef _DOUBLE
				sscanf_s(buffer,"%lf %lf",&pvalue,&vvalue);
#else
				sscanf_s(buffer,"%f %f",&pvalue,&vvalue);
#endif
				_pvtp_.push_back(pvalue);    //PVT压力(atm)
				_pvtv_.push_back(vvalue);    //PVT气体积系数
				++dcount;
			}
		}
		fclose(file);
	}
	return dcount;
}
void shiwen_connector::write_property(const string & fname,int savetype=0){
//	vector<REAL> tmp=vector<REAL>(_nnn_,SHIWEN_NAV);
	vector<REAL> tmp=vector<REAL>(_nnn_,_ncp_);
	vector<REAL>tmp2=vector<REAL>(_nnn_,0.0);
	if(savetype==0){ // by grids
		for(int i=0;i<_nn_;++i){
			tmp[i]=(REAL)_gridconnectormap_.cn[i];
			tmp2[i]=_oilabundance_[i];
		}
	}
	else{ // by nodes
		vector<int> indx=vector<int>();
		indx.reserve(8);
		int nindx=0;
		int k1=k2topz(_bottomlayer_)+1;
		int k2=k2topz(_toplayer_);
		for(int k=k1;k>=k2;--k){
			for(int j=0;j<_nny_;++j){
				for(int i=0;i<_nnx_;++i){
					nindx=k*_nnx_*_nny_+j*_nnx_+i;
					k2gindx(i,j,k,indx);
					REAL s=0.0;
#ifdef _DEBUG
	CString item;
	CString output;
	item.Format(_T("Node=(%d,%d,%d)"),i,j,k);
	output=item;
#endif
					for(size_t ii=0;ii<indx.size();++ii){
						s+=_oilabundance_[indx[ii]];
#ifdef _DEBUG
	int i1,j1,k1;
	indx2ijk(indx[ii],i1,j1,k1);
	item.Format(_T("Grid=(%d,%d,%d),Index=%d"),i1,j1,k1,indx[ii]);
	output+=item;
#endif
						if(_gridconnectormap_.cn[indx[ii]]>=0)
							tmp[nindx]=(REAL)(_gridconnectormap_.cn[indx[ii]]+1); // start from 1
					}
#ifdef _DEBUG
	item.Format(_T("s=%.6f\n"),s/indx.size());
	output+=item;
	OutputDebugString(output);
#endif
					tmp2[nindx]=s/indx.size();
				}
			}
		}		
	}
#ifdef _DEBUG
	write_property("cz3_tmp2.dat",&tmp2[0]);
#endif
	addproperty(SHIWEN_CONNECTOR_NUMBER,tmp);
	if(_type_==TYPE_OIL)
		addproperty(SHIWEN_ABUNDANCE,tmp2);
	else
		addproperty(SHIWEN_ABUNDANCEGAS,tmp2);
	OnWriteShiWenData(fname,*_pd_);
}
void shiwen_connector::write_property(const string &fname,REAL *pp,int bValid){
	FILE *pfile;
	int nindx, i,j,k;
	fopen_s(&pfile,fname.c_str(),"w");
	for(nindx=0;nindx<_nnn_;++nindx){
	//	if(bValid && !isvalid(nindx)) continue;
		if(bValid && !isvalid(pp,nindx)) continue;
		nindx2ijk(nindx,i,j,k);
	//	fprintf(pfile,"(%d,%d,%d)	%f\n",i,j,k,pp[nindx]);
		fprintf(pfile,"%f\n",pp[nindx]);
	}
	fclose(pfile);
}
void shiwen_connector::write_connector(const string &fname,const string &minfname){
	FILE *pfile;
	INT i;
//	INT i,n;
//	size_t ic;
	fopen_s(&pfile,fname.c_str(),"w");
	for(i=0;i<_nn_; ++i){
		REAL cn;
		if(_gridconnectormap_.cn[i]>=0)
			cn=(REAL)_gridconnectormap_.cn[i]+1;
		else
			cn=_ncp_;
		fprintf(pfile,"%f%s",cn,i<_nn_-1?"\n":"");
	}
	fclose(pfile);
/*	fopen_s(&pfile,minfname.c_str(),"w");
	fprintf(pfile,"ConnectorCount %d\n",_connector_.size());
	ic=0;
	for(vector<PCONNECTORGRAPH>::iterator it=_connector_.begin();it!=_connector_.end();++it,++ic){
		fprintf(pfile,"Connector %d(grids=%d)\n",ic,(n=(*it)->grids().size()));
		for(i=0;i<n;++i){
#ifdef _DEBUG
			int indx=(*it)->grids()[i]->indx,nindx=indx2nindx(indx);
			int ii,jj,kk;
			indx2ijk(indx,ii,jj,kk);
//			fprintf(pfile,"%d(%d, %d, %d) %f %f%s",indx,ii,jj,kk,_pd_->proValues[_zcoorpos_*_nnn_+nindx],
//				_pd_->proValues[_porositypos_*_nnn_+nindx], i<n-1||(i==n-1 && ic<_connector_.size()-1)?"\n":"");
			fprintf(pfile,"%d(%d, %d, %d) %f %f%s",indx,ii,jj,kk, _zcoor_[nindx],
				_porosity_[nindx], i<n-1||(i==n-1 && ic<_connector_.size()-1)?"\n":"");
#else
			fprintf(pfile,"%d%s",it->grids()[i].indx,i<n-1||(i==n-1 && ic<_connector_.size()-1)?"\n":"");
#endif
		}
	}
	fclose(pfile);*/
}
void shiwen_connector::write_oilabundance(const string &fname,const string &minfname){
	FILE *pfile;
	//size_t i,n, isl,j,pc;	
	int i;	
	fopen_s(&pfile,fname.c_str(),"w");
	i=0;
	for(vector<REAL>::iterator it=_oilabundance_.begin(); it!=_oilabundance_.end(); ++it,++i)
		fprintf(pfile,"%e%s",*it,i<_nn_-1?"\n":"");
	fclose(pfile);
/*
	fopen_s(&pfile,minfname.c_str(),"w");
	fprintf(pfile,"LinesCount %d\n",_pathnumberall_);

	isl=0;
	for(vector<PCONNECTORGRAPH>::iterator it=_connector_.begin();it!=_connector_.end();++it){
		for(i=0;i<(pc=(*it)->path().size());++i,++isl){
			fprintf(pfile,"Streamline %d(grids=%d)\n",isl,(n=(*it)->path()[i].size()));
			for(j=0;j<n;++j){
				int indx=(*it)->path()[i][j];
				int ii,jj,kk;
				indx2ijk(indx,ii,jj,kk);
				fprintf(pfile,"%f(%d, %d, %d)%s",_oilabundance_[indx],ii,jj,kk,j<n-1 ||(j==n-1) && isl<_pathnumberall_-1?"\n":"");
			}
		}
	}
	fclose(pfile);*/
}
void shiwen_connector::write_streamline(const string &fname){
	FILE *pfile;
	fopen_s(&pfile,fname.c_str(),"w");
	fprintf(pfile,"LinesCount %d\n",_pathnumberall_);
	size_t isl=1,i,j,pc,n;
	for(vector<PCONNECTORGRAPH>::iterator it=_connector_.begin();it!=_connector_.end();++it){
		for(i=0;i<(pc=(*it)->path().size());++i,++isl){
//			fprintf(pfile,"Streamline %d(grids=%d)\n",isl,(n=(*it)->path()[i].size()));
			fprintf(pfile,"Streamline %d %d\n",isl,(n=(*it)->path()[i].size()));
			for(j=0;j<n;++j){
//#ifdef _DEBUG
			/*int indx=it->path()[i][j];
			int ii,jj,kk;
			indx2ijk(indx,ii,jj,kk);
			fprintf(pfile,"%d(%d, %d, %d)%s",indx,ii,jj,kk,j<n-1 ||(j==n-1) && isl<_pathnumberall_-1?"\n":"");*/
			fprintf(pfile,"%d%s",(*it)->path()[i][j],j<n-1 ||(j==n-1) && isl<=_pathnumberall_-1?"\n":"");
/*#else
			fprintf(pfile,"%d%s",it->path()[i][j],j<n-1 ||(j==n-1) && isl<_pathnumberall_-1?"\n":"");
#endif*/
			}
			if(isl<_pathnumberall_)
				fprintf(pfile,"\n");
		}
	}
	fclose(pfile);
}
void shiwen_connector::write_streamline(const string &fname,const string &fname1){
	FILE *pfile,*pfile1;
	fopen_s(&pfile,fname.c_str(),"w");
	fopen_s(&pfile1,fname1.c_str(),"w");
#ifndef _DEBUG
	fprintf(pfile,"LinesCount %d\n",_pathnumberall_);
	fprintf(pfile1,"LinesCount %d\n",_pathnumberall_);
#endif
	size_t isl=1,i,j,pc,n;
	for(vector<PCONNECTORGRAPH>::iterator it=_connector_.begin();it!=_connector_.end();++it){
		for(i=0;i<(pc=(*it)->path().size());++i,++isl){
			n=(*it)->path()[i].size();
//			fprintf(pfile,"Streamline %d(grids=%d)\n",isl,(n=(*it)->path()[i].size()));
#ifndef _DEBUG
			fprintf(pfile,"Streamline %d %d\n",isl,n);
			fprintf(pfile1,"Streamline %d %d\n",isl,n);
#endif
			for(j=0;j<n;++j){
//#ifdef _DEBUG
			/*int indx=it->path()[i][j];
			int ii,jj,kk;
			indx2ijk(indx,ii,jj,kk);
			fprintf(pfile,"%d(%d, %d, %d)%s",indx,ii,jj,kk,j<n-1 ||(j==n-1) && isl<_pathnumberall_-1?"\n":"");*/
				int ii,jj,kk,indx;
				fprintf(pfile,"%d%s",(indx=(*it)->path()[i][j]),j<n-1 ||(j==n-1) && isl<=_pathnumberall_-1?"\n":"");
				indx2ijk(indx,ii,jj,kk);
			 	int nindx=ii+jj*_nx_+k2topz(kk)*_nx_*_ny_;
				fprintf(pfile1,"%d%s",nindx,j<n-1 ||(j==n-1) && isl<=_pathnumberall_-1?"\n":"");
/*#else
			fprintf(pfile,"%d%s",it->path()[i][j],j<n-1 ||(j==n-1) && isl<_pathnumberall_-1?"\n":"");
#endif*/
			}
#ifndef _DEBUG
			if(isl<_pathnumberall_){
				fprintf(pfile,"\n");
				fprintf(pfile1,"\n");
			}
#endif
		}
	}
	fclose(pfile);
	fclose(pfile1);
}

INT shiwen_connector::k2topz(INT k){
	int nl=0;
	for(int i=0;i<_pd_->ShiWenFileHead.nLayerNum;++i){
		nl+=(_pd_->vecEachLayerNumber[i]-1);
		if(nl>k){
			return k+i;
		}
	}
	return -1;
}
void shiwen_connector::k2gindx(int i,int j,int k,vector<int>& indx){
	indx.clear();
	int nl=0;
	int kk;
	for(kk=0;kk<_pd_->ShiWenFileHead.nLayerNum;++kk){
		nl+=_pd_->vecEachLayerNumber[kk];
		if(nl>k)
			break; // kk store large layernumber
	}
	int kg=k-kk; // layer of grid with current layer surface as top 
	int index0;
	if(k==nl-_pd_->vecEachLayerNumber[kk]) // top of ith large layer
	{
		if( i==0 && j==0){
			index0=kg*_nx_*_ny_;
			indx.push_back(index0);
		}
		if( i==0 && j==_ny_){
			index0=(j-1)*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
		}
		if( i==_nx_ && j==0){
			index0=i-1+kg*_nx_*_ny_;
			indx.push_back(index0);
		}
		if(i==_nx_ && j==_ny_){
			index0=i-1+(j-1)*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
		}
		if(i==0 && j>0 && j<_ny_){
			index0=i+(j-1)*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+j*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
		}
		if(i==_nx_ && j>0 && j<_ny_){
			index0=i-1+(j-1)*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
			index0=i-1+j*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
		}
		if(i>0 && i<_nx_ && j==0){
			index0=i-1+kg*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+kg*_nx_*_ny_;
			indx.push_back(index0);
		}	
		if(i>0 && i<_nx_ && j==_ny_){
			index0=i-1+(j-1)*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+(j-1)*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
		}
		if(i>0 && i<_nx_ && j>0 && j<_ny_){
			index0=i-1+(j-1)*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+(j-1)*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
			index0=i-1+j*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+j*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
		}
	}
	if(k==nl-1){
		if( i==0 && j==0){
			index0=(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
		}
		if( i==0 && j==_ny_){
			index0=(j-1)*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
		}
		if( i==_nx_ && j==0){
			index0=i-1+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
		}
		if(i==_nx_ && j==_ny_){
			index0=i-1+(j-1)*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
		}
		if(i==0 && j>0 && j<_ny_){
			index0=i+(j-1)*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+j*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
		}
		if(i==_nx_ && j>0 && j<_ny_){
			index0=i-1+(j-1)*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i-1+j*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
		}
		if(i>0 && i<_nx_ && j==0){
			index0=i-1+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
		}	
		if(i>0 && i<_nx_ && j==_ny_){
			index0=i-1+(j-1)*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+(j-1)*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
		}
		if(i>0 && i<_nx_ && j>0 && j<_ny_){
			index0=i-1+(j-1)*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+(j-1)*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i-1+j*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+j*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
		}
	}
	if(k>nl-_pd_->vecEachLayerNumber[kk]  && k<nl-1){
		if( i==0 && j==0){
			index0=(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=kg*_nx_*_ny_;
			indx.push_back(index0);
		}
		if( i==0 && j==_ny_){
			index0=(j-1)*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=(j-1)*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
		}
		if( i==_nx_ && j==0){
			index0=i-1+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i-1+kg*_nx_*_ny_;
			indx.push_back(index0);
		}
		if(i==_nx_ && j==_ny_){
			index0=i-1+(j-1)*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i-1+(j-1)*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
		}
		if(i==0 && j>0 && j<_ny_){
			index0=i+(j-1)*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+j*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+(j-1)*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+j*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
		}
		if(i==_nx_ && j>0 && j<_ny_){
			index0=i-1+(j-1)*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i-1+j*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i-1+(j-1)*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
			index0=i-1+j*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
		}
		if(i>0 && i<_nx_ && j==0){
			index0=i-1+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i-1+kg*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+kg*_nx_*_ny_;
			indx.push_back(index0);
		}	
		if(i>0 && i<_nx_ && j==_ny_){
			index0=i-1+(j-1)*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+(j-1)*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i-1+(j-1)*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+(j-1)*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
		}
		if(i>0 && i<_nx_ && j>0 && j<_ny_){
			index0=i-1+(j-1)*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+(j-1)*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i-1+j*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+j*_nx_+(kg-1)*_nx_*_ny_;
			indx.push_back(index0);
			index0=i-1+(j-1)*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+(j-1)*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
			index0=i-1+j*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
			index0=i+j*_nx_+kg*_nx_*_ny_;
			indx.push_back(index0);
		}
	}
}
int shiwen_connector::isvalid(PCONNECTORGRID pcgd){
	if( _zcoor_[pcgd->vertex.nindx111]!=SHIWEN_NAV && _zcoor_[pcgd->vertex.nindx211]!=SHIWEN_NAV
		&& _zcoor_[pcgd->vertex.nindx121]!=SHIWEN_NAV && _zcoor_[pcgd->vertex.nindx221]!=SHIWEN_NAV
		&& _zcoor_[pcgd->vertex.nindx112]!=SHIWEN_NAV && _zcoor_[pcgd->vertex.nindx212]!=SHIWEN_NAV
		&& _zcoor_[pcgd->vertex.nindx122]!=SHIWEN_NAV && _zcoor_[pcgd->vertex.nindx222]!=SHIWEN_NAV )
		return TRUE;
	else
		return FALSE;
}
void shiwen_connector::gridporosity(PCONNECTORGRID pcgd){
	pcgd->porosity=(_porosity_[pcgd->vertex.nindx111]+ _porosity_[pcgd->vertex.nindx211]+
		_porosity_[pcgd->vertex.nindx121]+ _porosity_[pcgd->vertex.nindx221]+
		_porosity_[pcgd->vertex.nindx112]+ _porosity_[pcgd->vertex.nindx212]+
		_porosity_[pcgd->vertex.nindx122]+ _porosity_[pcgd->vertex.nindx222])/8;
}
int shiwen_connector::isporosityvalid(PCONNECTORGRID pcgd){
	if(pcgd->porosity>=_porositylimit_/100.0)
		return TRUE;
	else
		return FALSE;
}
void shiwen_connector::gridoillostvolume(PCONNECTORGRID pcgd){
	if(_type_==TYPE_OIL)
		pcgd->lostvolume= (REAL)(pcgd->volume*pcgd->porosity*_residualoil_/100.0);		//残余油体积以有效孔隙体积为基准计算
	else{
		//溶解气体积以去除束缚水的体积为基准计算（与油计算方法不一致？）
		pcgd->lostvolume=pcgd->volume*_stoneden_*_adsorbedgas_;
		pcgd->lostvolume+=(REAL)(pcgd->volume*pcgd->porosity*(100.0-_constrainedwater_)/100.0*_residualoil_); 
	}
}
void shiwen_connector::gridsandstone(PCONNECTORGRID pcgd){
	pcgd->sandstone=(_sandstone_[pcgd->vertex.nindx111]+ _sandstone_[pcgd->vertex.nindx211]+
		_sandstone_[pcgd->vertex.nindx121]+ _sandstone_[pcgd->vertex.nindx221]+
		_sandstone_[pcgd->vertex.nindx112]+ _sandstone_[pcgd->vertex.nindx212]+
		_sandstone_[pcgd->vertex.nindx122]+ _sandstone_[pcgd->vertex.nindx222])/8;
}
void shiwen_connector::gridvalidporovolume(PCONNECTORGRID pcgd){
	if(_type_==TYPE_OIL)
		pcgd->validvolume=(REAL)(pcgd->volume*pcgd->porosity*( (100.0-_constrainedwater_)*pcgd->sandstone/1.0e4-_residualoil_/100.0) );
	else
	{

		REAL wf=(REAL)(_waterden_*_gravityacc_*pcgd->altitude*_pressurecoe_/1.e6*MP2ATM);
		REAL vc=p2vc(wf);
		pcgd->validvolume=(REAL)(pcgd->volume*pcgd->porosity*(100.0-_constrainedwater_)/100.0/vc-pcgd->lostvolume);
	}
}
void shiwen_connector::gridbuoyancy(PCONNECTORGRID pcgd)
{
	REAL h=pcgd->altitude;
	if(_type_==TYPE_OIL)
		pcgd->force_b=(REAL)((_waterden_-_oilden_)*_gravityacc_*h*_pressurecoe_);  //Pa 
	else{
		REAL wf=(REAL)(_waterden_*_gravityacc_*h*_pressurecoe_/1.e6*MP2ATM);
		REAL vc=p2vc(wf);
		pcgd->force_b=(REAL)((_waterden_-_oilden_/vc)*_gravityacc_*h*_pressurecoe_);  //Pa 
	}
}
void shiwen_connector::gridcapillary(PCONNECTORGRID pcgd){
	REAL r=(REAL)(EXP((REAL)((pcgd->porosity*100-17.228)/4.2707))/10000.0);   //um->cm
	pcgd->force_c=(REAL) (2.0*_tension_*NPM2KGPCM*COS(_anglecontact_*DEG2RAD)/r*KGPCM2TOMP*1.e6); // Pa
}

int shiwen_connector::gridproperty(PCONNECTORGRID pcgd,int i,int j,int k)
{
	int status=-1;
	pcgd->indx=i+j*_nx_+k*_nx_*_ny_; 
#ifdef _DEBUG
	if(pcgd->indx==199709)
		int a=1;
#endif
	pcgd->tz=k2topz(k);
	int nindx=i+j*_nnx_+pcgd->tz*_nnx_*_nny_;
	pcgd->vertex.nindx111=nindx;
	pcgd->vertex.nindx211=nindx+1;
	pcgd->vertex.nindx121=nindx+_nnx_;
	pcgd->vertex.nindx221=nindx+_nnx_+1;
	pcgd->vertex.nindx112=nindx+_nnx_*_nny_;
	pcgd->vertex.nindx212=nindx+_nnx_*_nny_+1;
	pcgd->vertex.nindx122=nindx+_nnx_*_nny_+_nnx_;
	pcgd->vertex.nindx222=nindx+_nnx_*_nny_+_nnx_+1;
	pcgd->accvolume=0.0;
	pcgd->isbottom=false;
	pcgd->isvisited=false;
	if(isvalid(pcgd)){
		status=0; // zcoor is valid
		gridporosity(pcgd);
		if(isporosityvalid(pcgd)){
			status=1; // zcoor and porosity are both valid
			gridsandstone(pcgd);
			pcgd->length=FABS((REAL)((_pd_->ShiWenFileHead.xMax-_pd_->ShiWenFileHead.xMin)/_nx_));		
			pcgd->width=FABS((REAL)((_pd_->ShiWenFileHead.yMax-_pd_->ShiWenFileHead.yMin)/_ny_));
			pcgd->height=FABS((REAL)((_zcoor_[pcgd->vertex.nindx112]-_zcoor_[pcgd->vertex.nindx111])
				+(_zcoor_[pcgd->vertex.nindx212]-_zcoor_[pcgd->vertex.nindx211])
				+(_zcoor_[pcgd->vertex.nindx122]-_zcoor_[pcgd->vertex.nindx121])
				+(_zcoor_[pcgd->vertex.nindx222]-_zcoor_[pcgd->vertex.nindx221]))/4);	
			INT ij=i+j*_nnx_;
			pcgd->altitude=FABS((REAL)((_zcoor_[pcgd->vertex.nindx112]+_zcoor_[pcgd->vertex.nindx111])
				+(_zcoor_[pcgd->vertex.nindx212]+_zcoor_[pcgd->vertex.nindx211])
				+(_zcoor_[pcgd->vertex.nindx122]+_zcoor_[pcgd->vertex.nindx121])
				+(_zcoor_[pcgd->vertex.nindx222]+_zcoor_[pcgd->vertex.nindx221]))/8)+
					(_surfaltitude_[ij]+_surfaltitude_[ij+1]+_surfaltitude_[ij+_nnx_]+_surfaltitude_[ij+_nnx_+1])/4;	
			pcgd->area=(REAL)((pcgd->length*pcgd->width+(pcgd->length+pcgd->width)*pcgd->height)*2/1e6); //(Km)^2
			pcgd->volume= (REAL)(pcgd->length*pcgd->width*pcgd->height/1e6);   //1e6 m^3 
			gridoillostvolume(pcgd);
			gridvalidporovolume(pcgd);
			gridbuoyancy(pcgd);
			gridcapillary(pcgd);
		}
		else
			return status;
	}
	else
		return status;
	return status;
}
int shiwen_connector::isvalid(INT nindx){
//	if(_pd_->proValues[_zcoorpos_*_nnn_+nindx] != SHIWEN_NAV)
	if(_zcoor_[nindx] != SHIWEN_NAV)
		return  TRUE;
	else
		return FALSE;
}
int shiwen_connector::isvalid(REAL *pp,INT nindx){
	if(pp[nindx] !=SHIWEN_NAV)
		return TRUE;
	else
		return FALSE;
}
int shiwen_connector::isporosityvalid(INT nindx){
//	if(_pd_->proValues[_porositypos_*_nnn_+nindx]>=_porositylimit_/100.0)
 	if(_porosity_[nindx]>=_porositylimit_/100.0)
		return TRUE;
	else
		return FALSE;
}
INTEGER shiwen_connector::ispropertyexisted(const string proname){
	vector<string>::iterator it=find(_pd_->proNames.begin(),_pd_->proNames.end(),proname);
	if(it == _pd_->proNames.end()) return FALSE;
	return it-_pd_->proNames.begin();
}
void shiwen_connector::addproperty(const string proname,vector<INT>& provalue){
	INTEGER pos;
	vector<REAL> tmp=vectori2r(provalue);
	if(!(pos=ispropertyexisted(proname))){//append new property
		_pd_->proNames .push_back(proname);
		_pd_->proValues.insert(_pd_->proValues.end(),tmp.begin(),tmp.end());
		_pd_->ShiWenFileHead.nProNum++;
	}
	else{//update
		copy(tmp.begin(),tmp.end(),_pd_->proValues.begin()+pos*_nnn_);
	}
}
void shiwen_connector::addproperty(const string proname,vector<REAL>& provalue){
	INTEGER pos;	
	vector<REAL>tmp=provalue;
	if(!(pos=ispropertyexisted(proname))){
		_pd_->proNames.push_back(proname);
		_pd_->proValues.insert(_pd_->proValues.end(),tmp.begin(),tmp.end());
		_pd_->ShiWenFileHead.nProNum++;
	}
	else{
		copy(tmp.begin(),tmp.end(),_pd_->proValues.begin()+pos*_nnn_);
	}
}
void shiwen_connector::rl2cstring(){
	/*
	int i,j,k=0;
	CString item;
	for(i=0;i<_pd_->ShiWenFileHead.nLayerNum;++i){
		for(j=0;j<_pd_->vecEachLayerNumber[i];++j){			
			item.Format(_T("%d  %s"),k,string2cstring(_pd_->vecEachLayerName[i]));
			_rlname_.push_back(item);
			++k;
		}
	}*/
	// 20190303
	CString item;
	_nz_=0;  // number of layers
	for(int i=0;i<_pd_->ShiWenFileHead.nLayerNum;++i){
		for(int j=0;j<_pd_->vecEachLayerNumber[i]-1;++j){
			item.Format(_T("%d-%d %s"),i,_nz_,string2cstring(_pd_->vecEachLayerName[i]));
			_rlname_.push_back(item);
			++_nz_;
		}
	}
}
void shiwen_connector::validgridnumallcount(){
	_validgridnumall_=0;
	for(int i=0;i<_nnn_;++i)
	{
//		if(_pd_->proValues[_zcoorpos_*_nnn_+i] != SHIWEN_NAV) 
		if(_zcoor_[i] != SHIWEN_NAV) 
			++_validgridnumall_;
	}
}		
INT shiwen_connector::nextvalidnindx(INT nindx){
	++nindx; // next node
	while(!isvalid(nindx))
		++nindx;
	if(nindx>_nnn_) nindx=-1;
	return nindx;
}
INT shiwen_connector::nextvalidnindx(REAL *pp,INT nindx){
	++nindx; // next node
	while(!isvalid(pp,nindx))
		++nindx;
	if(nindx>_nnn_) nindx=-1;
	return nindx;
}
int shiwen_connector::node2layer(){
	if(rlname().size()<3)
		return (int)rlname().size();
	_srclayer_--;
	_toplayer_--;
	_bottomlayer_--;
	if(_toplayer_<0)
		_toplayer_=0;
	if(_bottomlayer_<0)
		_bottomlayer_=0;
	if(_srclayer_<0)
		_srclayer_=0;
	return -100;
}
REAL shiwen_connector::gridarea(INT i, INT j, INT k){
		INT nindx=i+j*_nnx_+k*_nnx_*_nny_;
		REAL length=FABS((REAL)((_pd_->ShiWenFileHead.xMax-_pd_->ShiWenFileHead.xMin)/_nx_));		
		REAL width=FABS((REAL)((_pd_->ShiWenFileHead.yMax-_pd_->ShiWenFileHead.yMin)/_ny_));
//		REAL height=FABS((REAL)(_pd_->proValues[_zcoorpos_*_nnn_+nindx+_nnx_*_nny_]-_pd_->proValues[_zcoorpos_*_nnn_+nindx]));			
		REAL height=FABS((REAL)(_zcoor_[nindx+_nnx_*_nny_]-_zcoor_[nindx]));			

		return (REAL)((length*width+(length+width)*height)*2/1e6); //(Km)^2
}
REAL shiwen_connector::gridvolume(INT i, INT j, INT k){
//	INT  nindx=i+j*_nnx_+k*_nnx_*_nny_;
	INT nindx=i+j*_nnx_+k2topz(k)*_nnx_*_nny_;
	REAL length=FABS((REAL)((_pd_->ShiWenFileHead.xMax-_pd_->ShiWenFileHead.xMin)/_nx_));		
	REAL width=FABS((REAL)((_pd_->ShiWenFileHead.yMax-_pd_->ShiWenFileHead.yMin)/_ny_));
//	REAL height=FABS((REAL)(_pd_->proValues[_zcoorpos_*_nnn_+nindx+_nnx_*_nny_]-this->_pd_->proValues[_zcoorpos_*_nnn_+nindx]));			
	REAL height=FABS((REAL)(_zcoor_[nindx+_nnx_*_nny_]-_zcoor_[nindx]));			

	return (REAL)(length*width*height/1e6);   //1e6 m^3 
}
REAL shiwen_connector::porovolume(INT i, INT j, INT k){
	INT nindx=i+j*_nnx_+k*_nnx_*_nny_;
//	return (REAL)( gridvolume(i,j,k)*this->_pd_->proValues[_porositypos_*_nnn_+nindx]/1.e2);
	return (REAL)( gridvolume(i,j,k)*_porosity_[nindx]);
}

REAL shiwen_connector::oillostvolume(INT i, INT j, INT k){
	if(_type_==TYPE_OIL)
		return (REAL)(porovolume(i,j,k)*_residualoil_/100.0);		//残余油体积以有效孔隙体积为基准计算
	else{
		REAL vol=gridvolume(i,j,k);
		REAL lost=vol*_stoneden_*_adsorbedgas_;
		INT nindx=i+j*_nnx_+k*_nnx_*_nny_;
		lost+=(REAL)(porovolume(i,j,k)*(100.0-_constrainedwater_)/100.0*_residualoil_); //溶解气体积以去除束缚水的体积为基准计算（与油计算方法不一致？）
		return lost;
	}
}
REAL shiwen_connector::validporovolume(INT i, INT j, INT k){
	INT nindx=i+j*_nnx_+k*_nnx_*_nny_;
	if(_type_==TYPE_OIL)
		return (REAL)(porovolume(i,j,k)*( (100.0-_constrainedwater_)*_sandstone_[nindx]/1.0e4-_residualoil_/100.0) );
	else
	{
		INT ij=i+j*_nx_;
//		REAL h=(REAL)(0.5*(this->pd()->proValues[_zcoorpos_*_nnn_+nindx]+this->pd()->proValues[_zcoorpos_*_nnn_+nindx+_nnx_*_nny_])+this->_surfaltitude_[ij]);
		REAL h=(REAL)(0.5*(_zcoor_[nindx]+_zcoor_[nindx+_nnx_*_nny_])+this->_surfaltitude_[ij]);
		REAL wf=(REAL)(_waterden_*_gravityacc_*h*_pressurecoe_/1.e6*MP2ATM);
		REAL vc=p2vc(wf);
		return (REAL)(porovolume(i,j,k)*(100.0-_constrainedwater_)/100.0/vc-oillostvolume(i,j,k));
	}
}
void shiwen_connector::indx2ijk(INT indx, INT &i, INT &j,INT &k){
		k=indx/((this->_pd_->ShiWenFileHead.nXNodeNum-1)*(this->_pd_->ShiWenFileHead.nYNodeNum-1));
		INT ij=indx-k*(this->_pd_->ShiWenFileHead.nXNodeNum-1)*(this->_pd_->ShiWenFileHead.nYNodeNum-1);
		j=ij/(this->_pd_->ShiWenFileHead.nXNodeNum-1);
		i=ij-j*(this->_pd_->ShiWenFileHead.nXNodeNum-1);
}

void shiwen_connector::nindx2ijk(INT nindx, INT &i,INT &j, INT &k){
		k=nindx/(_nnx_*_nny_);
		INT ij=nindx-k*_nnx_*_nny_;
		j=ij/_nnx_;
		i=ij-j*_nnx_;
}
INT shiwen_connector::indx2nindx(int indx){
	int i,j,k;
	indx2ijk(indx,i,j,k);
	return i+j*_nnx_+k*_nnx_*_nny_;
}

REAL shiwen_connector::throatsize(INT i, INT j, INT k){
	INT  nindx=i+j*_nnx_+k*_nnx_*_nny_;
//	return (REAL)(EXP((REAL)((this->_pd_->proValues[_porositypos_*_nnn_+nindx]-17.228)/4.2707))/10000.0);   //um->cm
	return (REAL)(EXP((REAL)((_porosity_[nindx]*100-17.228)/4.2707))/10000.0);   //um->cm

}
REAL shiwen_connector::capillary_pressure(INT i,INT j, INT k){
	return (REAL) (2.0*_tension_*NPM2KGPCM*COS(_anglecontact_*DEG2RAD)/throatsize(i,j,k)*KGPCM2TOMP*1.e6);
}
REAL shiwen_connector::calc_weight(INT indx1, int indx2){
	/*  weight is negtive of driving force, represent aibilty of oil from indx1 to indx2
	*/
	REAL bf1,bf2,cp1,cp2;
		// driving force
/*	calc_f(indx1,bf1,cp1);
	calc_f(indx2,bf2,cp2);*/
	// oil will go to higher grid with lower boyrant force
	// or go to grid with larger porosity with lower capillary force 
	// however weight is negtive of driving force
	// we assume that oil goes from indx1 to indx2, then bouyancy driving force is bf2-bf1(<0, will to indx2)
	// capillary driving force is cp2-cp1(<0, will to indx2)
	bf1=_cgds_[_gridconnectormap_.pos[indx1]].force_b;
	cp1=_cgds_[_gridconnectormap_.pos[indx1]].force_c;
	bf2=_cgds_[_gridconnectormap_.pos[indx2]].force_b;
	cp2=_cgds_[_gridconnectormap_.pos[indx2]].force_c;
	return -(bf1-bf2+cp1-cp2);  
}
void shiwen_connector::calc_f(INT indx,REAL &bf,REAL &cp){
	// static boyrant force is drving force,
	//projection into coordinate faces, including surface
	INT i,j,k,ij,nindx; 
	k=indx/(_nx_*_ny_);
	ij=indx-k*_nx_*_ny_;
	j=ij/_nx_;
	i=ij-j*_nx_;
	nindx=i+j*_nnx_+k*_nnx_*_nny_;/**/

	// buoyant force
//	REAL h=(REAL)(0.5*(this->pd()->proValues[_zcoorpos_*_nnn_+nindx]+this->pd()->proValues[_zcoorpos_*_nnn_+nindx+_nnx_*_nny_])+this->_surfaltitude_[ij]);
//	REAL h=(REAL)(0.5*(_zcoor_[nindx]+_zcoor_[nindx+_nnx_*_nny_])+this->_surfaltitude_[ij]);
	REAL h=_cgds_[_gridconnectormap_.pos[indx]].altitude;
	if(_type_==TYPE_OIL)
		bf=(REAL)((_waterden_-_oilden_)*h*_pressurecoe_);  //Pa 
	else{
		REAL wf=(REAL)(_waterden_*_gravityacc_*h*_pressurecoe_/1.e6*MP2ATM);
		REAL vc=p2vc(wf);
		bf=(REAL)((_waterden_-_oilden_/vc)*_gravityacc_*h*_pressurecoe_);  //Pa 
	}
	// capillary pressure
	k=_cgds_[_gridconnectormap_.pos[indx]].tz;
	cp=capillary_pressure(i,j,k);  //Pa
}
REAL shiwen_connector::p2vc(REAL p){
	vector<REAL>::iterator it=lower_bound(_pvtp_.begin(),_pvtp_.end(),p);
	if(it==_pvtp_.end()) return _pvtp_.back();
	if(it==_pvtp_.begin()) return _pvtp_.front();
	INTEGER i=it-_pvtp_.begin();
	return _pvtv_[i-1]+(_pvtv_[i]-_pvtv_[i-1])/(_pvtp_[i-1]-_pvtp_[i])*(_pvtp_[i-1]-p);
}
REAL shiwen_connector::vc2p(REAL vc){
	vector<REAL>::iterator it=lower_bound(_pvtv_.begin(),_pvtv_.end(),vc,greater<REAL>());
	if(it==_pvtv_.begin()) return _pvtv_.front();
	if(it==_pvtv_.end()) return _pvtv_.back();
	INTEGER i=it-_pvtv_.begin();
	return _pvtp_[i-1]+(_pvtp_[i]-_pvtp_[i-1])/(_pvtv_[i-1]-_pvtv_[i])*(_pvtv_[i-1]-vc);
}
void shiwen_connector::edge(PCONNECTORGRID pcgd,INT preindx){
	// weight is negtive, means that oil goes from indx to preindx
	REAL weight=calc_weight(pcgd->indx,preindx); 

	// path of preindx
	PWEIGHTEDEDGE pwg= new WEIGHTEDEDGE(pcgd->indx,-weight);
	INT pos=_gridconnectormap_.pos[preindx];
	vector<WEIGHTEDEDGE>::iterator it1=lower_bound(_cgds_[pos].edge.begin(),_cgds_[pos].edge.end(),*pwg,AscendingSortByWeight());	
	_cgds_[pos].edge.insert(it1,*pwg);
	delete pwg;

	// path of indx
	pwg=new WEIGHTEDEDGE(preindx,weight);
	pos=_gridconnectormap_.pos[pcgd->indx];
	vector<WEIGHTEDEDGE>::iterator it2=lower_bound(_cgds_[pos].edge.begin(),_cgds_[pos].edge.end(),*pwg,AscendingSortByWeight());
    _cgds_[pos].edge.insert(it2,*pwg);
	delete pwg;
}
void shiwen_connector::connectorgraphcreate(PCONNECTORGRID pcgd){
	this->_connector_.push_back(new connectorgraph(pcgd));
}
void shiwen_connector::connectorgraphpop()
{
	_connector_.back()->clear();
	_connector_.pop_back();
}
void shiwen_connector::connectorgraphtraversebfs(){
	INT posnc=-1,pos,i,ii,status;
	INT gridsize;
	vector<bool> isvisited((gridsize=(INT)_cgds_.size()),false);
	queue<INT> q;
	while(++posnc<gridsize){
		if(isvisited[posnc]) continue;
		q.push(posnc);  // start of graph
		isvisited[posnc]=true;
		status=0;
		while(! q.empty()){
			i=q.front();	// take first element
			// here are operations on the first element
			if(_cgds_[i].isbottom) status=1;
			if(i==posnc)
				this->connectorgraphcreate(&_cgds_[i]);
			else
				this->_connector_.back()->append(&_cgds_[i]);
			q.pop();			// remove first element
			for(vector<WEIGHTEDEDGE>::iterator it=_cgds_[i].edge.begin();it!=_cgds_[i].edge.end();++it)
				if(!isvisited[pos=_gridconnectormap_.pos[it->indx]]){
					q.push(pos);
					isvisited[pos]=true;
			}
		}
		if(!status) 
			connectorgraphpop();  // not valid because on source
		else{ // create connector-grid map table
			sort(_connector_.back()->grids().begin(),_connector_.back()->grids().end(),DescendingSortByIndex());
			ii=0;
			for(vector<PCONNECTORGRID>::iterator it=_connector_.back()->grids().begin();
				it !=_connector_.back()->grids().end();++it){
				_gridconnectormap_.cn[(*it)->indx]=(int)_connector_.size()-1;
				_gridconnectormap_.pos[(*it)->indx]=ii++;
			}
		}
	}
}
int shiwen_connector::checkpropertypos(int format){
	int status=0;
	if(_zcoorpos_== ZCOOR_POS) 
		status+=ZCOOR_POS;
	else if(_zcoorpos_<_pd_->ShiWenFileHead.nProNum)
	{
		_zcoor_.clear();
		for(int i=0;i<_nnn_;++i)
			_zcoor_.push_back(_pd_->proValues[_zcoorpos_*_nnn_+i]);
		validgridnumall();
	}		
	if(_porositypos_==POROSITY_POS)
		status+=POROSITY_POS;
	else if(_porositypos_<_pd_->ShiWenFileHead.nProNum)
	{
		_porosity_.clear();
#ifdef _DEBUG
		int a=1;
#endif
		for(int i=0;i<_nnn_;++i){
			_porosity_.push_back(_pd_->proValues[_porositypos_*_nnn_+i]);
			if(_maxporosity_<_porosity_[i])
				_maxporosity_=_porosity_[i];
#ifdef _DEBUG
			if(_porosity_[i]!=SHIWEN_NAV)
				a=2;			
#endif
		}
#ifdef _DEBUG
		CString item;
		item.Format(_T("a=%d"),a);
		OutputDebugString(item);
#endif
	}
	if(format!=0){ // %
		_maxporosity_=0.0;
		for(int i=0;i<_nnn_;++i){
			_porosity_[i]/=100.0;
			if(_maxporosity_<_porosity_[i])
				_maxporosity_=_porosity_[i];
		}
	}
	if(_oilgenpos_ == OILGEN_POS)
		status+=OILGEN_POS;
	else if(_oilgenpos_<_pd_->ShiWenFileHead.nProNum)
	{
		_oilgen_.clear();
		for(int i=0;i<_nnn_;++i)
			_oilgen_.push_back(_pd_->proValues[_oilgenpos_*_nnn_+i]);
	}
	return status;
}
int shiwen_connector::find_connector(){
	int status=-1,statustmp=-1; // no valid grid
	INT indx,preindx;
	INT i,j,k,valid=FALSE,NONE=TRUE,gridsize=0;
	CONNECTORGRID cgd;
	// 初始化
	flushconnector();

	//int cdcn=0;	
	/*	connection direction of current node
				  	0, no connector
				    1, backward, 2, right, 4, down; 
				  	3, backward & right
				 	5, backward & down
				 	6, right & download
				 	7, backward, right, down */
	
	//linear search, O(N), loop over meshgrid
	for(k=_bottomlayer_;k>=_toplayer_;--k){ // bottom layer number starting from 0
		if(k<_bottomlayer_ && NONE==TRUE) {return 0;} // no valid contact,return
		for(j=_ny_-1;j>=0;--j){
			for(i=_nx_-1;i>=0;--i){
				if((statustmp=gridproperty(&cgd,i,j,k))>0){
					status=1; // at least one connector
					if(k==_bottomlayer_) {cgd.isbottom=true;NONE=FALSE;}
					// put current grid into buffer and record its position
					_cgds_.push_back(cgd);	
					_gridconnectormap_.pos[indx=cgd.indx]=gridsize++;
					if(i<_nx_-1) // backward:
						if (_gridconnectormap_.pos[preindx=indx+1]>-1)
							this->edge(&cgd,preindx);
					if(j<_ny_-1) // right
						if (_gridconnectormap_.pos[preindx=indx+_nx_]>-1)
							this->edge(&cgd,preindx);
					if(k<_bottomlayer_) //down
						if (_gridconnectormap_.pos[preindx=indx+_nx_*_ny_]>-1)
							this->edge(&cgd,preindx);


			/*		if(i<_nx_-1) if (_gridconnectormap_.pos[indx+1]>-1)	cdcn+=1;				// backward:
					if(j<_ny_-1) if (_gridconnectormap_.pos[indx+_nx_]>-1)	cdcn+=2;				// right
					if(k<_bottomlayer_) if (_gridconnectormap_.pos[indx+_nx_*_ny_]>-1) cdcn+=4;	// down		

					switch(cdcn){
					case 0: //new connection			
						break;
					// 1,2,4 added into exited _connector_
					case 1:
						preindx=indx+1;
						this->edge(&cgd,preindx);
						break;
					case 2:
						preindx=indx+_nx_;
						this->edge(&cgd,preindx);
						break;
					case 4:
						preindx=indx+_nx_*_ny_;
						this->edge(&cgd,preindx);
						break;
					// 3,5,6,7 connect two or three exited connectors	
					case 3:	
						preindx=indx+1;
						this->edge(&cgd,preindx);

						preindx=indx+_nx_;
						this->edge(&cgd,preindx);
						break;
					case 5:
						preindx=indx+1;
						this->edge(&cgd,preindx);

						preindx=indx+_nx_*_ny_;
						this->edge(&cgd,preindx);
						break;
					case 6:
						preindx=indx+_nx_;
						this->edge(&cgd,preindx);

						preindx=indx+_nx_*_ny_;
						this->edge(&cgd,preindx);
						break;
					case 7:
						preindx=indx+1;
						this->edge(&cgd,preindx);

						preindx=indx+_nx_;
						this->edge(&cgd,preindx);
 
						preindx=indx+_nx_*_ny_;
						this->edge(&cgd,preindx);
						break;
					default:
						// some error message can be printed
						break; 
					}*/
				}
				else if(statustmp==0){ // at least one's zcoor is valid
					if(_minimaxporosity_<cgd.porosity)
						_minimaxporosity_=cgd.porosity;
					status=0;
				}
			}
		}
	}
	// put grids into connector buffer
	connectorgraphtraversebfs();
	if(_connector_.size()>0)
		status=(int)_connector_.size();
	return status;
}

REAL shiwen_connector::oildischargeamount(INT i, INT j, INT k){
/*	INT  nindx=i+j*_nnx_+k*_nnx_*_nny_;
//	REAL oilgen=this->pd()->proValues[_oilgenpos_*_nnn_+nindx];
	REAL oilgentmp=_oilgen_[nindx];
	oilgentmp*=(REAL)(_oildischargecoe_/100.0*this->gridarea(i,j,k));  
	if(oilgentmp<_oilgenmin_) oilgentmp=(REAL) 0.0;
	if(_type_==TYPE_OIL)
		return (REAL)(oilgentmp/_oilden_*10.0);	// 1e6 m^3
	else
	{
		INT ij=i+j*_nx_;
//		REAL h=(REAL)(0.5*(this->pd()->proValues[_zcoorpos_*_nnn_+nindx]+this->pd()->proValues[_zcoorpos_*_nnn_+nindx+_nnx_*_nny_])+this->_surfaltitude_[ij]);
		REAL h=(REAL)(0.5*(_zcoor_[nindx]+_zcoor_[nindx+_nnx_*_nny_])+this->_surfaltitude_[ij]);
		REAL wf=(REAL)(_waterden_*_gravityacc_*h*_pressurecoe_/1.e6*MP2ATM);
		REAL vc=p2vc(wf);
		return (REAL)(oilgentmp/vc/1.e2);
	}*/
	int indx=i+j*_nx_+k*_nx_*_ny_; 
	int tz=k2topz(k);
	int nindx=i+j*_nnx_+tz*_nnx_*_nny_;
	int nindx111=nindx;
	int nindx211=nindx+1;
	int nindx121=nindx+_nnx_;
	int nindx221=nindx+_nnx_+1;
	int nindx112=nindx+_nnx_*_nny_;
	int nindx212=nindx+_nnx_*_nny_+1;
	int nindx122=nindx+_nnx_*_nny_+_nnx_;
	int nindx222=nindx+_nnx_*_nny_+_nnx_+1;

	REAL length=FABS((REAL)((_pd_->ShiWenFileHead.xMax-_pd_->ShiWenFileHead.xMin)/_nx_));		
	REAL width=FABS((REAL)((_pd_->ShiWenFileHead.yMax-_pd_->ShiWenFileHead.yMin)/_ny_));
	if(_zcoor_[nindx112]==SHIWEN_NAV || _zcoor_[nindx111] ==SHIWEN_NAV ||
		_zcoor_[nindx212]==SHIWEN_NAV ||_zcoor_[nindx211]==SHIWEN_NAV ||
		_zcoor_[nindx122]==SHIWEN_NAV ||_zcoor_[nindx121]==SHIWEN_NAV ||
		_zcoor_[nindx222]==SHIWEN_NAV ||_zcoor_[nindx221]==SHIWEN_NAV)
		return 0.0;
	if(_oilgen_[nindx112]==SHIWEN_NAV || _oilgen_[nindx111] ==SHIWEN_NAV ||
		_oilgen_[nindx212]==SHIWEN_NAV ||_oilgen_[nindx211]==SHIWEN_NAV ||
		_oilgen_[nindx122]==SHIWEN_NAV ||_oilgen_[nindx121]==SHIWEN_NAV ||
		_oilgen_[nindx222]==SHIWEN_NAV ||_oilgen_[nindx221]==SHIWEN_NAV)
		return 0.0;
	REAL height=FABS((REAL)((_zcoor_[nindx112]-_zcoor_[nindx111])
		+(_zcoor_[nindx212]-_zcoor_[nindx211])
		+(_zcoor_[nindx122]-_zcoor_[nindx121])
		+(_zcoor_[nindx222]-_zcoor_[nindx221]))/4);	
	REAL area=(REAL)((length* width+( length+ width)* height)*2/1e6); //(Km)^2
	REAL oilgentmp=(_oilgen_[nindx111]+ _oilgen_[nindx211]+
		_oilgen_[nindx121]+ _oilgen_[nindx221]+
		_oilgen_[nindx112]+ _oilgen_[nindx212]+
		_oilgen_[nindx122]+ _oilgen_[nindx222])/8;
#ifdef _DEBUG
		if(maxoilgen<oilgentmp)
			maxoilgen=oilgentmp;
#endif
		oilgentmp*=(REAL)(_oildischargecoe_/100.0*area);  
		if(oilgentmp<_oilgenmin_) oilgentmp=(REAL) 0.0;
		if(_type_==TYPE_OIL)
			return (REAL)(oilgentmp/_oilden_*10.0);	// 1e6 m^3
		else{
			INT ij=i+j*_nnx_;
			REAL h=(REAL)FABS((REAL)((_zcoor_[nindx112]+_zcoor_[nindx111])
				+(_zcoor_[nindx212]+_zcoor_[nindx211])
				+(_zcoor_[nindx122]+_zcoor_[nindx121])
				+(_zcoor_[nindx222]+_zcoor_[nindx221]))/8)+
				(_surfaltitude_[ij]+_surfaltitude_[ij+1]+_surfaltitude_[ij+_nnx_]+_surfaltitude_[ij+_nnx_+1])/4;	
			REAL wf=(REAL)(_waterden_*_gravityacc_*h*_pressurecoe_/1.e6*MP2ATM);
			REAL vc=p2vc(wf);
			return (REAL)(oilgentmp/vc/1.e2);    
		}
}
REAL shiwen_connector::oildischargeamount(INT indx){
	INT i,j,k;
	indx2ijk(indx,i,j,k);
	return oildischargeamount(i,j,k);
}

int shiwen_connector::oillost(INT indx,REAL &lostamount,REAL &srcamount){
	PCONNECTORGRID pcgd=this->_connector_[_gridconnectormap_.cn[indx]]->grids()[_gridconnectormap_.pos[indx]];
	REAL tmp=pcgd->lostvolume;
	if(tmp<=srcamount){ //source is enough
		pcgd->lostvolume=0;		// only lost once
		lostamount+=tmp;
		srcamount-=tmp;
		if(srcamount<=(REAL)0.0) return FALSE;
	}
	else{
		pcgd->lostvolume-=srcamount;
		lostamount+=srcamount;
		srcamount=(REAL) 0.0;
		return FALSE;
	}
	return TRUE;
}
int shiwen_connector::findminweight(INT indx){
	INT cn;
	int newindx=-1;
	PCONNECTORGRID pcgd0=_connector_[(cn=_gridconnectormap_.cn[indx])]->grids()[_gridconnectormap_.pos[indx]],pcgd;
	if(pcgd0->validvolume<=(REAL)0.0){ // current grid is full, which means oil can goes along positive weight
		for(vector<WEIGHTEDEDGE>::iterator it=pcgd0->edge.begin();it!=pcgd0->edge.end();++it){
			pcgd=_connector_[cn]->grids()[_gridconnectormap_.pos[it->indx]];
		//	if(!pcgd->isvisited){
			//	pcgd->isvisited=true;
				if(pcgd->validvolume>(REAL)0.0){// weight is least, not visited and have valid volume
					newindx=it->indx;
					break;
				}
			//}
		}
	}
	else{ // current grid is not full
		for(vector<WEIGHTEDEDGE>::iterator it=pcgd0->edge.begin();it!=pcgd0->edge.end();++it){
			pcgd=_connector_[cn]->grids()[_gridconnectormap_.pos[it->indx]];
		//	if(!pcgd->isvisited){
			//	pcgd->isvisited=true;
				if(pcgd->validvolume>(REAL)0.0){
					if(it->weight<0.0) // only negtive weight is valid
						newindx=it->indx;
					break; 
				}
		//	}
		}
	}
	// only zero valid volume (i.e. full) can block advance path, so two different advance path could share nodes
	// for new advance path, the visted nodes are agained not visited
/*	for(vector<WEIGHTEDEDGE>::iterator it=pcgd0->edge.begin();it!=pcgd0->edge.end();++it){
		_connector_[cn]->grids()[_gridconnectormap_.pos[it->indx]]->isvisited=false;
	}*/
	return newindx;//-1 indicate, no valid next grid, may be bound, or all neighbours are visited
}
int shiwen_connector::findmaxweight(INT indx){
	INT cn;
	size_t i;
	PCONNECTORGRID pcgd=_connector_[(cn=_gridconnectormap_.cn[indx])]->grids()[_gridconnectormap_.pos[indx]];	
	for(i=pcgd->edge.size();i>0;--i){
		indx=pcgd->edge[i-1].indx; //indx changed
		if(this->_connector_[cn]->grids()[_gridconnectormap_.pos[indx]]->validvolume>(REAL)0.0)
			return indx;
	}
	return -1;
}
int shiwen_connector::advancepath(PCONNECTORGRAPH pgraph,INT &indx,REAL & srcamount,REAL &lostamount,vector<INT>&path){
	// assume source is enough, if not, this advance path is invalid
	int status=TRUE;
	INT nextindx;
	lostamount=(REAL)0.0; // for every path, its initial loss is zero
	path.clear();

	pgraph->grids()[_gridconnectormap_.pos[indx]]->isvisited=true;  // label grid as visited
	path.push_back(indx);
	if(oillost(indx,lostamount,srcamount)){ // source enought
		for(nextindx=findminweight(indx); nextindx>=0;nextindx=findminweight(indx)){//next valid grid
			indx=nextindx;
			path.push_back(indx);
			if(!oillost(indx,lostamount,srcamount)){
				status=FALSE;
				break;
			}
		}
		if(nextindx<0){ // no valid path,fill
			status=TRUE;
		}
	}
	else{
		status=FALSE;
	}
	// only zero valid volume (i.e. full) can block advance path, so two different advance path could share nodes
	// for new advance path, the visted nodes are agained not visited
/*	for(vector<PCONNECTORGRID>::iterator it=pgraph->grids().begin();it !=pgraph->grids().end(); ++it)
		(*it)->isvisited=false;*/
	return status;
} 
int shiwen_connector::fillgrid(INT indx,REAL & srcamount){
	INT cn;
	REAL tmp;
	PCONNECTORGRID pcgd=_connector_[(cn=_gridconnectormap_.cn[indx])]->grids()[_gridconnectormap_.pos[indx]];
	if((tmp=srcamount-pcgd->validvolume)<=(REAL)0.0){ // source is not enough
		pcgd->validvolume-=srcamount;
		pcgd->accvolume+=srcamount;
		this->_connector_[cn]->setvalidvolume(this->_connector_[cn]->validvolume()-srcamount);
		srcamount=(REAL)0.0;
		return FALSE;
	}
	else{ // source enough
		this->_connector_[cn]->setvalidvolume(this->_connector_[cn]->validvolume()-pcgd->validvolume);
		pcgd->accvolume+=pcgd->validvolume;
		pcgd->validvolume=(REAL)0.0;
		srcamount=tmp;
		return TRUE;
	}
}
void shiwen_connector::migrationinconnector(PCONNECTORGRAPH pgraph){
	INT indx;
	REAL srcamountall, srcamount,lostamountall,lostamount;
	size_t n=pgraph->grids().size();
	vector<INT> path;	//path buffer
	srcamountall=(REAL)0.0;
	lostamount=(REAL)0.0;
	lostamountall=(REAL)0.0;

	for(size_t ii=0; ii<n && pgraph->grids()[ii]->isbottom;++ii){ // source contact node loop
		indx=pgraph->grids()[ii]->indx;
		srcamount=oildischargeamount(indx+_nx_*_ny_);  // source amount
		srcamountall+=srcamount;
		
		if(pgraph->path().empty() ||(!pgraph->path().empty() && !pgraph->path().back().empty()))pgraph->createpath(); // only create new path when need
		int ap;
#ifdef _DEBUG
		if(indx==199709)
			int a=1;
#endif 
		while((ap=advancepath(pgraph,indx,srcamount,lostamount,path))){// find accumulation point
			lostamountall+=lostamount;	
			pgraph->appendpath(path);
			if(!fillgrid(indx,srcamount)) break; // no enough source		
			if((indx=findmaxweight(indx))<0) break; // a potential boundary, at least sourroudning nodes is full
		}
		if(ap==FALSE){
			lostamountall+=lostamount;	
			pgraph->appendpath(path);
		}
	}
#ifdef _DEBUG
		CString item;
		item.Format(_T("Max Oil Gen Intensity %.6f\n"), maxoilgen);
		OutputDebugString(item);
#endif
	if(pgraph->path().back().empty()) pgraph->poppath();
	// statistics
	pgraph->setsrcamount(srcamountall);
	pgraph->setoillost(lostamountall);
	// oil abundance: 1e4t /1e6m^3
/*	REAL vol;
	for(vector<PCONNECTORGRID>::iterator it=pgraph->grids().begin();it!=pgraph->grids().end();++it){
		if((vol=(*it)->volume)<=REAL_EPSILON)
			vol=REAL_EPSILON;
		if(_type_==TYPE_OIL)
			_oilabundance_[(*it)->indx]=((*it)->accvolume*_oilden_)/vol*KGPM3TO4TP6M3;  //1e^4t /1e6m^3

		else
			_oilabundance_[(*it)->indx]=(REAL)(((*it)->accvolume/1e2)/vol);  // 1e8m^3/1e6m^3
	}*/
	REAL area;
	for(vector<PCONNECTORGRID>::iterator it=pgraph->grids().begin();it!=pgraph->grids().end();++it){
		if((area=(*it)->area)<=REAL_EPSILON)
			area=REAL_EPSILON;
		if(_type_==TYPE_OIL)
			_oilabundance_[(*it)->indx]=((*it)->accvolume*_oilden_)/area*KG6TO4T;  //1e^4t /1e6m^3

		else
			_oilabundance_[(*it)->indx]=(REAL)(((*it)->accvolume/1e2)/area);  // 1e8m^3/1e6m^3
	}

}
int shiwen_connector::migration(){
	_pathnumberall_=0;
	for(vector<PCONNECTORGRAPH>::iterator it=_connector_.begin();it!=_connector_.end();++it)
	{
		migrationinconnector(*it);
		_pathnumberall_+=(*it)->path().size();
	}
	return (int)_pathnumberall_;
}

#ifdef _DEBUG
/*int shiwen_connector::test_connector(int mn)
{
	pPeroModData pd=_pd_;
	INT i,j,k,nx=pd->ShiWenFileHead.nXNodeNum,ny=pd->ShiWenFileHead.nYNodeNum,nz=pd->ShiWenFileHead.nZNodeNum,nn=nx*ny*nz,n;
	REAL * pp,* po,* po1,* pz=&(pd->proValues[0]);
	int status=0;
	int xmin ,xmax ,ymin ,ymax ,zmin ,zmax ;
	int xc ,yc ;			
	int sx ,sy ,sz ;// source contact node
	int kk=0;
	errno_t  err;
	
	xmin=70<nx-1?70:0,xmax=80<nx-1?80:nx-1,ymin=100<ny-1?100:0,ymax=120<ny-1?120:ny-1,zmin=10<nz-1?10:0,zmax=18<nz-1?18:nz-1;
	xc=(xmin+xmax)/2,yc=(ymin+ymax)/2;			
	sx=121<nx-5?121:nx-5,sy=123<ny-5?123:ny-5,sz=18<nz-4?18:nz-4;// source contact node
	_oilgenpos_=_oilgenpos_>_pd_->ShiWenFileHead.nProNum?_pd_->ShiWenFileHead.nProNum-1:_oilgenpos_;
	_porositypos_=_porositypos_>_pd_->ShiWenFileHead.nProNum?_pd_->ShiWenFileHead.nProNum-1:_porositypos_;

	switch(mn){
	case 0:
		break;
	case 1:
		pp=&(pd->proValues[_porositypos_*nn]);
		po=&(pd->proValues[_oilgenpos_*nn]);
		// initialization
		for(i=0;i<nn;++i)
			pp[i]=SHIWEN_NAV;
		srand((int)time(0));
		for(i=xmin;i<=xmax;++i)
			for(j=ymin;j<=ymax;++j){
				for(k=zmin;k<=zmax;++k){
						pp[i+j*nx+k*nx*ny]=RANDOM(8.0,12.0);
				}
			}
		break;
	case 2:
		pp=&(pd->proValues[_porositypos_*nn]);
		po=&(pd->proValues[_oilgenpos_*nn]);
	// initialization
	for(i=0;i<nn;++i)
		pp[i]=SHIWEN_NAV;
		kk=0;
		for(i=xmin;i<=xc;++i)
			for(j=ymin;j<=yc;++j)
				for(k=zmin;k<=zmax;++k){
					++kk;
					pp[i+j*nx+k*nx*ny]=(REAL)8;
				}
		kk=0;
		for(i=xc+2;i<=xmax;++i)
			for(j=yc+2;j<=ymax;++j)
				for(k=zmin;k<=zmax;++k){
					++kk;
					pp[i+j*nx+k*nx*ny]=(REAL)8;
				}
		break;
	case 3:
		pp=&(pd->proValues[_porositypos_*nn]);
		po=&(pd->proValues[_oilgenpos_*nn]);
	// initialization
	for(i=0;i<nn;++i)
		pp[i]=SHIWEN_NAV;
		pp[xmin+ymin*nx+zmax*nx*ny]=(REAL) 8; 
		pp[xmin+(ymin+1)*nx+zmax*nx*ny]=(REAL) 8; 
		pp[xmin+(ymin+2)*nx+zmax*nx*ny]=(REAL) 8; 

		pp[xmin+1+(ymin+2)*nx+zmax*nx*ny]=(REAL) 8; 
		pp[xmin+1+(ymin+3)*nx+zmax*nx*ny]=(REAL) 8;

		pp[xmin+2+ymin*nx+zmax*nx*ny]=(REAL) 8; 
		pp[xmin+2+(ymin+1)*nx+zmax*nx*ny]=(REAL) 8; 
		pp[xmin+2+(ymin+2)*nx+zmax*nx*ny]=(REAL) 8; 
		pp[xmin+2+(ymin+4)*nx+zmax*nx*ny]=(REAL) 8; 

		pp[xmin+3+(ymin+2)*nx+zmax*nx*ny]=(REAL) 8; 
		pp[xmin+3+(ymin+3)*nx+zmax*nx*ny]=(REAL) 8; 
		pp[xmin+3+(ymin+4)*nx+zmax*nx*ny]=(REAL) 8; 
		break;
	case 4:
		pp=&(pd->proValues[_porositypos_*nn]);
		po=&(pd->proValues[_oilgenpos_*nn]);
	// initialization
	for(i=0;i<nn;++i)
		pp[i]=SHIWEN_NAV;
			//(0,j,0)
		for(i=xmin,j=ymin,k=zmin;j<=ymax;++j)
			pp[i+j*nx+k*nx*ny]=(REAL) 8; 
		for(i=xmin,j=ymax,k=zmin;i<=xmax;++i)
			pp[i+j*nx+k*nx*ny]=(REAL) 8; 
		for(i=xmax,j=ymax,k=zmin;k<=zmax;++k)
			pp[i+j*nx+k*nx*ny]=(REAL) 8; 
		for(i=xmax,j=ymax,k=zmax;j>=ymin;--j)
			pp[i+j*nx+k*nx*ny]=(REAL) 8; 
		for(i=xmax,j=ymin,k=zmax;i>=xmin;--i)
			pp[i+j*nx+k*nx*ny]=(REAL) 8;
		break;
	case 5:
		pp=&(pd->proValues[_porositypos_*nn]);
		po=&(pd->proValues[_oilgenpos_*nn]);
	// initialization
	for(i=0;i<nn;++i)
		pp[i]=SHIWEN_NAV;
		pp[ sx+nx*sy+nx*ny*sz			]=RANDOM(16.5,17.5); 
		pp[ sx+nx*sy+nx*ny*(sz-1)		]=RANDOM(15.5,16.5); 

		pz[ sx+nx*sy+nx*ny*(sz-1)		]=(REAL)((pz[sx+nx*sy+nx*ny*sz]+pz[sx+nx*sy+nx*ny*(sz-2)])*0.5);

		pp[ sx+nx*sy+nx*ny*(sz-2)		]=RANDOM(11.5,14.0);  //center
		pp[ sx-2+nx*sy+nx*ny*(sz-2)		]=RANDOM(11.5,14.0); 
		pp[ sx-1+nx*sy+nx*ny*(sz-2)		]=RANDOM(11.5,14.0); 
		pp[ sx+1+nx*sy+nx*ny*(sz-2)		]=RANDOM(11.5,14.0); 
		pp[ sx+2+nx*sy+nx*ny*(sz-2)		]=RANDOM(11.5,14.0); 
		pp[ sx+nx*(sy-2)+nx*ny*(sz-2)	]=RANDOM(11.5,14.0); 
		pp[ sx+nx*(sy-1)+nx*ny*(sz-2)	]=RANDOM(11.5,14.0); 
		pp[ sx+nx*(sy+1)+nx*ny*(sz-2)	]=RANDOM(11.5,14.0); 
		pp[ sx+nx*(sy+2)+nx*ny*(sz-2)	]=RANDOM(11.5,14.0);

		pp[ sx+nx*sy+nx*ny*(sz-3)		]=RANDOM(8.5,9.5); 
		pz[ sx+nx*sy+nx*ny*(sz-3)		]=(REAL)((pz[ sx+nx*sy+nx*ny*(sz-2)]+pz[ sx+nx*sy+nx*ny*(sz-4)])*.5);				
		pp[ sx+nx*sy+nx*ny*(sz-4)		]=RANDOM(7.5,8.5); 

 		po[sx+nx*sy+nx*ny*(sz+1)]=1700;
		break;
	case 6:
		_porositypos_=9;
		_oilgenpos_=10;
		pp=&(pd->proValues[_porositypos_*nn]);
		po=&(pd->proValues[_oilgenpos_*nn]);
		if(!porosity_mul){
			for(i=0;i<nn;++i)
				pp[i]*=(REAL)100.0;
			porosity_mul=1;
		}
		//sandstone
		for(i=0;i<_nnn_;++i)
			_sandstone_[i]=(REAL)((pd->proValues[1*_nnn_+i]+pd->proValues[2*_nnn_+i])/(pd->proValues[1*_nnn_+i]+pd->proValues[2*_nnn_+i]+pd->proValues[3*_nnn_+i]+
			+pd->proValues[4*_nnn_+i]+pd->proValues[5*_nnn_+i]+pd->proValues[6*_nnn_+i]+pd->proValues[7*_nnn_+i]+pd->proValues[8*_nnn_+i])*100.0);
		//read oilgenensity
		FILE *pfile;
		err=fopen_s(&pfile,"ML_38_HQ_19_6N_1_15layer_eclipse - oilgenintensity.txt","r");
		if(err) 
			return mn; // error happened
		po1=po;
		n=(nx-1)*(ny-1)*(nz-1);
		for(i=0;i<n;++i)
			fscanf_s(pfile,"%f",po1++);
		fclose(pfile);
		break;
	default:
		break;
	}	
	return status;
}*/
size_t mem;
size_t GetMemMB(){
	HANDLE handle = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS pmc;
	GetProcessMemoryInfo(handle, &pmc, sizeof(pmc));
	return pmc.WorkingSetSize/1024/1024;
	//printf("%d\r\n",pmc.WorkingSetSize); 
}
#endif
