#pragma once

// 头文件
#include "ShiWen_iPeraMod_DataInterface.h"
#include<algorithm>
#include<queue>
#include<stack>
#include<functional>
#include<stdio.h>

//数据类型
typedef iPeraModData * pPeroModData;
#ifdef X86
typedef int INTEGER;
#else
typedef _int64 INTEGER;
#endif
typedef int	INT;
#ifdef _DOUBLE
typedef double REAL;
#define REALCSTRINGFORMAT "%lf"
#define REAL_MAX DBL_MAX
#undef REAL_EPSILON
#define REAL_EPSILON DBL_EPSILON
#else
typedef float REAL ;
#define REALCSTRINGFORMAT "%f"
#define REAL_MAX FLT_MAX
#undef REAL_ESPSION
#undef REAL_EPSILON
#define REAL_EPSILON FLT_EPSILON
#endif
#undef REAL_MIN
#define REAL_MIN -REAL_MAX

//宏定义
#define SUPERTRUE 2
#define TRUE 1
#define FALSE 0
#define TYPE_NONE	-1
#define TYPE_OIL	0
#define TYPE_GAS	1
#define POROSITY_POS -1		// porosity property position
#define OILGEN_POS	-2		// oil or gas generating amount position
#define ZCOOR_POS -4		// z coordinate positions
#define ABUNDANCE_OIL "oilabundance"
#define ABUNDANCE_GAS "gasabundance"
//默认参数值与单位转换常数
const REAL SHIWEN_NAV=(REAL)-99999.0;							// 可考虑作为整个工程的全局变量，可以一致的表示无效值
const REAL SHIWEN_NCP=-1;										// 非连通网格属性值
const REAL MATH_PI=(REAL)3.1415926535897932384626433832795;		// pi
const REAL GRAVITY_ACC=(REAL) 9.81;								// m/s^2
const REAL ANGLE_CONTACT=(REAL)0.00;							// 0
const REAL WATER_DENSITY=(REAL)1000.0;							// kg/m^3
const REAL KG6TO4T=(REAL)0.1;									// Rho(Kg/m^3)  =   Rho(1e4 t/1e6 m^3)
const REAL KGPM3TO4TP6M3=(REAL)0.1;								// Rho(Kg/m^3)  =   Rho(1e4 t/1e6 m^3)
const REAL DEG2RAD=(REAL)(MATH_PI/180.0);						// angle(rad)	=	angle(deg)*DEG2RAD
const REAL RAD2DEG=(REAL)(180.0/MATH_PI);						// angle(deg)	=	angle(rad)*RAD2DEG
const REAL MP2ATM=(REAL)9.86923;								// p(atm)		=	p(MPa)*MP2ATM
const REAL N2KG=(REAL)(1.0/GRAVITY_ACC);						// f(kg)		=	f(N)*N2KG
const REAL NPM2KGPCM=(REAL) (N2KG/100.0);						// sigma(N/m)	=	sigma(Kg/cm)
const REAL KGPCM2TOMP=(REAL)(0.01*GRAVITY_ACC);					// p(MPa)		=	p(kg.cm^-2)*KGPCM2TOMP	
const REAL SHIWEN_POROSITY=(REAL)7.0;							// valid porocity, %
const REAL SHIWEN_TENSION=(REAL)0.0175;						// oil water interfacial tension, N/m
const REAL SHIWEN_OIL_DENSITY=(REAL)900.0;						// oil density, kg/m^3
const REAL SHIWEN_SANDSTONE_PC=(REAL)70.0;						// sandstone percent content, %
const REAL SHIWEN_PRESSURE_COE=(REAL)1.2;						// pressure coefficience, ps=ph*SHIWEN_PRESSURE_COE
const REAL SHIWEN_VALID_POROCITY_COE=(REAL)0.5;					// valid porocity, Vv=Vp*SHIWEN_VALID_POROCITY_COE 
const REAL SHIWEN_OIL_DISCHARGE_COE=(REAL)50.0;					// oil discharge coefficient, %
const REAL SHIWEN_OIL_GEN_INTENSITY_MIN=(REAL)0.0;				// minimum oil generating intensity
const REAL SHIWEN_SURFACE_ALTITUDE=200.0;						// surface altilude
const REAL SHIWEN_CONTRAINEDWATER=(REAL)30.0;					// constrained water saturation
const REAL SHIWEN_RESIDUALOIL=(REAL)5.0;						// residual oil saturation
const REAL SHIWEN_GASDEN=(REAL)1.036;							// gas density
const REAL SHIWEN_STONEDEN=(REAL)2.5;							// stone specific weight, t/m^3
const REAL SHIWEN_ADSORBEDGAS=(REAL)0.01;						// adsorbed gas, m^3/t
const REAL SHIWEN_DISSOLVEDGAS=(REAL)1.0;						// gas amount in water, m^3/m^3
const REAL SHIWEN_GAS_DISCHARGE_COE=(REAL)90.0;					// gas discharge coefficient, %
const REAL SHIWEN_GAS_GEN_INTENSITY_MIN=(REAL)0.0;				// minimum gas generating intensity
const REAL SHIWEN_PERSSURE_ERROR=(REAL)0.0;						// allowed pressure error (atm）
const string SHIWEN_CONNECTOR_NUMBER="connectorNumber";
const string SHIWEN_ABUNDANCE="abundanceOil";
const string SHIWEN_ABUNDANCEGAS="abundanceGas";


// 通用函数声明
REAL EXP(REAL);
REAL COS(REAL);
REAL FABS(REAL);
REAL RANDOM();				//产生随机数：0-1
REAL RANDOM(REAL, REAL);	//指定上下界的随机数
vector<REAL> vectori2r(vector<INT>&);	//向量数据类型转换
vector<INT> vectorr2i(REAL *,int n);	
CString string2cstring(string str);		//C++ string to CString
CString int2cstring(int n);

// 数据类型:
//     连通体由图描述, 图的顶点为网格, 网格相邻构成边, 边的权值为驱动力的相反数
typedef struct weigthededge{//边
	INT indx;
	REAL weight;
	weigthededge(INT x, REAL y){indx=x;weight=y;}
} WEIGHTEDEDGE, *PWEIGHTEDEDGE;

struct AscendingSortByWeight{//边按权（油）升序
// here is the function that will be called by std:sort/low_bound/upper_bound
	bool operator()(const WEIGHTEDEDGE &lhs, const WEIGHTEDEDGE &rhs){
		if( lhs.weight < rhs.weight) return true;
		if( rhs.weight < lhs.weight) return false;
		return lhs.weight < rhs.weight;
	}
};
typedef struct corner{
	INT nindx111,nindx211,nindx121,nindx221,nindx112,nindx212,nindx122,nindx222; // 8-corner points
}CORNER;
typedef struct connectorgrid{//网格(顶点)
	INT indx;
	CORNER vertex;
	INT tz;									// top z index
	REAL width;
	REAL length;
	REAL height;
	REAL area;
	REAL volume;							// 1e6 m^3
	REAL altitude;
	REAL porosity;
	REAL sandstone;
	REAL validvolume;						// valid volume:excluding non-sandstone, non-porosity, constrained water, residual oil
	REAL lostvolume;						// lost amount
	REAL accvolume;							// accumulation volume
	bool isbottom;							// contact source
	bool isvisited;							// indicate weather it is visited or not
	REAL force_b;							// buoyancy force, Pa
	REAL force_c;							// capillary pressure,Pa
	vector<WEIGHTEDEDGE> edge;				// edge within one connector

	REAL chargeoilamount;					// source property
} CONNECTORGRID, *PCONNECTORGRID;

struct DescendingSortByIndex{//网格（顶点）按编号降序
// Here is the function that will be called by std:sort/low_bound/upper_bound
    bool operator()( const PCONNECTORGRID &lhs, const PCONNECTORGRID &rhs ) {
        if( lhs->indx > rhs->indx ) return true;
        if( rhs->indx > lhs->indx ) return false;
        return lhs->indx > rhs->indx;
    }
};

typedef class connectorgraph{//连通体

	REAL _volume_;			//连通体体积
	REAL _validvolume_;		//有效体积
	REAL _oillost_;			//油损失量（体积单位）
	REAL _srcamount_;		//烃源油量（体积单位）
	vector<PCONNECTORGRID> _grids_;		//网格向量
	vector<vector<INT>> _path_;			//运聚路径

public:	
	connectorgraph();//构造函数
	connectorgraph(PCONNECTORGRID pcgd);

	REAL volume(){return _volume_;}				//属性返回值
	REAL validvolume(){return _validvolume_;}
	REAL oillost(){return _oillost_;}
	REAL srcamount(){return _srcamount_;}
	vector<PCONNECTORGRID> &grids(){return _grids_;}		//引用，返回向量本身
	vector<vector<INT>>& path(){return _path_;}

	void setvalidvolume(REAL vv){_validvolume_=vv;}		//设置属性值
	void setoillost(REAL ol){_oillost_=ol;}
	void setsrcamount(REAL sa){_srcamount_=sa;}

	void append(PCONNECTORGRID pcgd);	//添加网格
	void createpath();					//运聚路径创建
	void appendpath(vector<INT>& path);	//运聚路径添加
	void poppath();						//删除最后路径

	void clear();	//清空内存
}CONNECTORGRAPH,*PCONNECTORGRAPH;

class shiwen_connector{ //网格整体类
	int _type_;					// 运聚类型：TYPE_OIL/TYPE_GAS
	pPeroModData _pd_;			// 网格数据指针
	
    int _zcoorpos_,_oilgenpos_,_porositypos_; // 属性位置
	INT _nx_,_ny_,_nz_,_nn_,_nnx_,_nny_,_nnz_,_nnn_;	//网络维度
	INT _validgridnumall_;		// 有效网格总数
	INT _toplayer_;				// 储层顶层
	INT _bottomlayer_;			// 储层底层
	INT _srclayer_;				// 烃源层
	REAL _ncp_;					// 非连通网格属性值
	REAL _pressurecoe_;			// 压力系数
	REAL _porositylimit_;		// 有效孔隙度(%)
	REAL _constrainedwater_;	// 束缚水饱和度(%)
	REAL _waterden_;			// 水密度(Kg/m^3)
	REAL _oilden_;				// 油密度(Kg/m^3)、				气密度(Kg/m^3)
	REAL _gravityacc_;			// 重力加速度（m/s^2)
	REAL _tension_;				// 表面张力（N/m）
	REAL _anglecontact_;		// 接触角
	REAL _oildischargecoe_;		// 排油系数（%）,				生气系数（%）
	REAL _oilgenmin_;			// 最低有效生油量,				最低生气量
	REAL _pressureallowederr_;	//								压力允许误差 (atm)
	REAL _residualoil_;			// 残余油饱和度 (%), 			溶解气（m^/m^3)
	REAL _stoneden_;			//								岩石比重 t/m^3,仅气用
	REAL _adsorbedgas_;			//								吸附气 m^3/t，仅气用
	size_t _pathnumberall_;		// 路径总数
	vector<REAL> _surfaltitude_;		// 地表海拔
	vector<REAL> _pvtp_;				// pvt曲线
	vector<REAL> _pvtv_;				// pvt 曲线
	vector<CString> _rlname_;			// 储层名字（编号+地层名）
	vector<REAL> _porosity_;
	REAL _maxporosity_;
	REAL _minimaxporosity_;
	int _porosityformat_;
	vector<REAL> _oilgen_;				
	vector<REAL> _zcoor_;
	vector<REAL> _sandstone_;			// 砂岩百分含量(%)
	int _ssstatus_;
	vector<CONNECTORGRID> _cgds_;		// 网格缓存

	struct{
		vector<INT> cn;
		vector<INT> pos;
	}_gridconnectormap_;					// 网格-连通体编号映射表

	vector<PCONNECTORGRAPH> _connector_;	// 连通体
	vector<REAL> _oilabundance_;			// 资源丰度(10^7Kg/10^6m^3,万吨/平方公里*米)		或	(10^8m^3/10^6m^3，亿方/平方公里*米)
//	string _abundancepropertyname_;			// 资源丰度名称

public:
	shiwen_connector(int t=TYPE_OIL);				// 构造函数
	~shiwen_connector();			// 析构函数
	int initialize();				// 初始化:0初始化错误，1初始化正确
	void reinit();					// 二次初始化,仅初始化油气区别的参数
	void flushconnector();			// 重新初始化连通体内存
	void clear();					// 清空内存 
	int checkpropertypos(int format);			// 检测默认属性位置

	//	与MFC接口，用户自定义主函数可参，读写属性数据文件
	int read_sproperty(const string &,vector<REAL>&);
	int read_property(const string &);
	int read_pvt(const string &);
	void write_property(const string &,int);
	void write_property(const string &,REAL *,int bValid=FALSE);
	void write_connector(const string &,const string &minfname);
	void write_oilabundance(const string &,const string &minfname);
	void write_streamline(const string &);
	void write_streamline(const string &,const string &);

	int type(){return _type_;}			// 返回属性值
	int zcoorpos(){return _zcoorpos_;}
	int oilgenpos(){return _oilgenpos_;}
	int porositypos(){return _porositypos_;}
	pPeroModData pd(){return _pd_;}
	INT nx(){return _nx_;}
	INT ny(){return _ny_;}
	INT nz(){return _nz_;}
	INT nn(){return _nn_;}
	INT nnx(){return _nnx_;}
	INT nny(){return _nny_;}
	INT nnz(){return _nnz_;}
	INT nnn(){return  _nnn_;}	
	INT validgridnumall(){return _validgridnumall_;}
	INT toplayer(){return _toplayer_;}
	INT bottomlayer(){return _bottomlayer_;}
	INT srclayer(){return _srclayer_;}
	REAL ncp(){return _ncp_;}
	REAL pressurecoe(){return _pressurecoe_;}
	REAL porositylimit(){return _porositylimit_;}
	REAL constrainedwater(){return _constrainedwater_;}
	REAL waterden(){return _waterden_;}
	REAL oilden(){return _oilden_;}
	REAL gravityacc(){return _gravityacc_;}
	REAL tension(){return _tension_;}
	REAL anglecontact(){return _anglecontact_;}
	REAL oildischargecoe(){return _oildischargecoe_;}
	REAL oilgenmin(){return _oilgenmin_;}
	REAL residualoil(){return _residualoil_;}
	REAL stoneden(){return _stoneden_;}
	REAL adsorbedgas(){return _adsorbedgas_;}
	REAL pressureallowederr(){return _pressureallowederr_;}
	REAL maxporosity(){return _maxporosity_;}
	REAL minimaxporosity(){return _minimaxporosity_;}
	int porosityformat(){return _porosityformat_;}
	size_t pathnumberall(){return _pathnumberall_;}
	vector<CString> & rlname(){return _rlname_;}
	vector<REAL> &surfaltitude(){return _surfaltitude_;}
	vector<REAL>& oilabundance(){return _oilabundance_;}
	vector<REAL>& sandstone(){return _sandstone_;}
	int ssstatus(){return _ssstatus_;}
	vector<REAL>& porosity(){return _porosity_;};
	vector<REAL>& oilgen(){return _oilgen_;};				
	vector<REAL>& zcoor(){return  _zcoor_;}
	vector<PCONNECTORGRAPH> &connector(){return this->_connector_;}

	void settype(int t){_type_=t;}							//设置储藏类型
	void setzcoorpos(int zp){_zcoorpos_=zp;}
	void setoilgenpos(int gp){_oilgenpos_=gp;}
	void setporositypos(int pp){_porositypos_=pp;}
	void settoplayer(INT tl){_toplayer_=tl;}				//设置属性
	void setbottomlayer(INT bl){_bottomlayer_=bl;}
	void setsrclayer(INT sl){_srclayer_=sl;}
	void setncp(REAL n){_ncp_=n;}
	void setpressurecoe(REAL opc){_pressurecoe_=opc;}
	void setporosity(REAL p){_porositylimit_=p;}
	void setconstrainedwater(REAL cw){_constrainedwater_=cw;}
	void setsandstone(REAL ss){vector<REAL>(_nnn_,ss).swap(_sandstone_);}	
	void setssstatus(){_ssstatus_=1;}
	void unsetsstatus(){_ssstatus_=0;}
	void setsurfaltitude(vector<REAL> &sa){_surfaltitude_=sa;}
	void setwaterden(REAL wd){_waterden_=wd;}
	void setoilden(REAL od){_oilden_=od;}
	void setgravityacc(REAL ga){_gravityacc_=ga;}
	void settension(REAL t){_tension_=t;}
	void setanglecontact(REAL ac){_anglecontact_=ac;}
	void setoildischargecoe(REAL odc){_oildischargecoe_=odc;}
	void setoilgenmin(REAL ogm){_oilgenmin_=ogm;}
	void setresidualoil(REAL ro){_residualoil_=ro;}
	void setstoneden(REAL sd){_stoneden_=sd;}
	void setadsorbedgas(REAL ag){_adsorbedgas_=ag;}
	void setpressureallowederr(REAL pae){_pressureallowederr_=pae;}
	void setporosityformat(int format){_porosityformat_=format;}

	INT k2topz(INT k);
	void k2gindx(int i,int j,int k,vector<int>& indx);
	int isvalid(PCONNECTORGRID pcgd);
	void gridporosity(PCONNECTORGRID pcgd);
	int isporosityvalid(PCONNECTORGRID pcgd);
	void gridoillostvolume(PCONNECTORGRID pcgd);
	void gridsandstone(PCONNECTORGRID pcgd);
	void gridvalidporovolume(PCONNECTORGRID pcgd);
	void gridbuoyancy(PCONNECTORGRID pcgd);
	void gridcapillary(PCONNECTORGRID pcgd);
	int gridproperty(PCONNECTORGRID pcgd,int i,int j,int k);

	int isvalid(INT nindx);									// 判断网格属性是否有效（根据ZCOOR属性）:0无效，1有效
	int isvalid(REAL *pp,INT nindx);
	int isporosityvalid(INT nindx);							// 检测孔隙度是否有效:0无效，1有效
	INTEGER  ispropertyexisted(const string proname);			//检测属性是否存在：0不存在，1存在
	void addproperty(const string proname,vector<INT> &);	//添加属性,根据属性内容重载，当前所有属性用实型
	void addproperty(const string proname,vector<REAL> &);
	void rl2cstring(); 										// 获得层位名称向量:_rlname_
	void validgridnumallcount();							// 统计有效网格数
	INT nextvalidnindx(INT nindx);							// 查找下一个有效网格
	INT nextvalidnindx(REAL *pp,INT nindx);
	int node2layer();

	REAL gridarea(INT i, INT j, INT k);			// 计算网格面积(Km)^2
	REAL gridvolume(INT i, INT j, INT k);		// 计算网格体积(Km)^3->1e6 m^3
	REAL porovolume(INT,INT,INT);				// 孔隙体积，考虑孔隙度,1e6 m^3,
	REAL oillostvolume(INT,INT,INT);			// 油损失体积，只考虑一次;		气损失量
	REAL validporovolume(INT,INT,INT);			// 有效体积，去除束缚水，油（气）损失量

	INT indx2nindx(INT indx);							// 网格编号转节点编号
	void indx2ijk(INT indx, INT &i, INT &j,INT &k);		// 网络线性编号转三维编号
	void nindx2ijk(INT nindx, INT &i,INT &j, INT &k);	// 节点线性编号转三维编号 

	REAL throatsize(INT i, INT j,INT k);				// 经验喉道半径， cm
	REAL capillary_pressure(INT i, INT j, INT k);		// 毛细管压力，P=2*sigma*cos(theta)/r, Pa
	REAL calc_weight(INT indx1, INT indx2);				// 驱动力权重，为驱动力负数，油气运聚选择最小权重路径
	void calc_f(INT indx,REAL &bf,REAL &cp);			// 计算驱动力，地层压力与毛细管力

	REAL vc2p(REAL vc);									// 体积系数转压力
	REAL p2vc(REAL p);									// 压力转体积系数

	// 连通体操作接口
	void edge(PCONNECTORGRID pcgd,INT preindx);		// 添加连通体路径
	void connectorgraphcreate(PCONNECTORGRID pcgd);				// 创建新的连通体
	void connectorgraphpop();									// 删除最后一个连通体
	void connectorgraphtraversebfs();							// 遍历与组装连通体
	int find_connector();										// 查找连通体,返回连通体个数 
	int oillost(INT indx,REAL &lostamount,REAL &srcamount);		// 网格吸附油量计算,返回值：0,烃源油量不足;1，烃源油充足 
	REAL oildischargeamount(INT i,INT j,INT k);					// 计算烃源油量（1e6 m^3) 
	REAL oildischargeamount(INT indx);				
	REAL gasdischargeamount(INT i,INT j,INT k);					// 计算烃源气量（1e6 m^3) 
	REAL gasdischargeamount(INT indx);				
	int findminweight(INT indx);								// 查找最小权重：仅负权重有效，假设油不能逆驱动力运移
	int findmaxweight(INT indx);								// 查找最大权重：聚集可以受迫进行，故假设可以逆驱动力聚集		
	int advancepath(PCONNECTORGRAPH pgraph,INT &indx,
		REAL & srcamount,REAL &lostamount,vector<INT> &path);	// 以当前节点为起点，查找下一个聚集点，返回值：0,烃源不足；烃源充足
	int fillgrid(INT,REAL &srcamount);							// 聚集点网格油填充,返回值：0，烃源不足;1，烃源充足
	void migrationinconnector(PCONNECTORGRAPH pgraph);			// 当前连通体运聚
	int migration();											// 运聚接口函数，调用migrationinconnector函数，返回路径数

#ifdef _DEBUG
/*  
	mn indicate model number:
  		mn=1, simple 11*21*9 cubic connector at the top-left-back corner
  		mn=2, two connectors, 6*11*9 cubic one in model 1 and another 4*9*9 one
		mn=3, one 2d connector with irregular shape at layer 18 

  					    0 1 2  3  4
  	               0    . . .
     	           1	    .  .
     	           2	. . .     .     
      	           3	    .  .  .

		mn=4, one connector with irregular shape within 11*21*9 cubic, 
			  first on the top, y-axis (70,j,10) --> x-direct (i,120,10)
			  then z-direct (80,120,k), --> negtive y-direct (80,j,18)
			  then negtive x-direct (i,100,18)
					   --------
				    		  /
				     	     /
				     /	    |
				    /       |
	               ----------
        mn=5, a 3d asterisk with soucre  at (121,123,19)
			            o    *
			            o  *                  
		            - - o - -
		              * o
	                *   o

		mn=6, a real model
*/
int test_connector(int mn=0);
#endif
};
#ifdef _DEBUG
#include<Psapi.h>
#pragma comment(lib,"Psapi.lib")
extern size_t mem;
size_t GetMemMB();
#endif