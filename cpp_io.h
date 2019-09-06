////////////////////////////////////////////////////////////////
// 版权所有 (c) 2017, 石文软件有限公司
// 保留所有权利
// 
// 文件名称：ShiWen_iPeraMod_DataInterface.h
// 摘    要：
//			盆模系统网格模型文件头
//
// 使用方法：		
// 
// 平　台：Windows 9x/NT/2k/XP/2003
// 编译器：Visual C++ .Net
// 
// 当前版本：1.0
// 作者    ：刘少兵
//
// 开始日期：2017-08-16  
//			
// 完成日期：
// 修改    ： 孙亮-RIPED
//
// 开始日期：2017-09-21
// 最后修改：2017-09-26
////////////////////////////////////////////////////////////////

#pragma once
#include <string>
#include <vector>
#include <map>
#include <iostream>  
#include <fstream>
using namespace std;
//typedef unsigned short UINT;
//add by Ainray, 20180822, for NAN
typedef std::numeric_limits<float> Info;
#define NAN (Info::quiet_NaN());
//石文数据体文件头： by 石文
struct ShiWenFileHead
{
	ShiWenFileHead()
	{
		nFileFlag = 0x8001;//盆模文件格式标记
		xMin = 0;
		xMax = 0;
		yMin = 0;
		yMax = 0;
		zMin = 0;
		zMax = 0;
		nXNodeNum = 0;
		nYNodeNum = 0;
		nZNodeNum = 0;
		nProNum = 0;
		nProNameLength = 120;//120个字节
		nLayerNum = 0; //初始化为0个层
	}
	int    nFileFlag;      //文件标识名
	double xMin;           //X轴最小值
	double xMax;           //X轴最大值
	double yMin;           //Y轴最小值
	double yMax;           //Y轴最大值
	double zMin;           //Z轴最小值
	double zMax;           //Z轴最大值
	int     nXNodeNum;     //X轴节点个数
	int     nYNodeNum;     //Y轴节点个数
	int     nZNodeNum;     //Z轴节点个数
	int     nProNum;       //属性名称个数
	int     nProNameLength;//属性名称固定长度
	int     nLayerNum;     //层的个数,文件紧接着存nLayerNum个 int型 的层序号
};

//注意节点的书写方式：顺序是Z-Y-X，也就是说最外圈是Z！！！！！！！！！！！！

//iPeramod盆地模拟数据接口 by 孙亮 from RIPED
struct iPeraModData
{
	//第一部分：网格数据，采用石文网格格式
	ShiWenFileHead ShiWenFileHead;
	//第二部分：属性数据，采用附加数据格式
	//map<string, vector<float> > pros; // 修改：float 换成 double
	vector<string> proNames; // string 存储属性名称
//	vector<vector<float> > proValues; //vector为属性值
	vector<float> proValues;
	//第三部分：附加信息
	vector<int> vecEachLayerNumber; //每个大层的小层层数，为了从同一个文件能识别出多个地层
	vector<string> vecEachLayerName; //每个大层的层名 by xuwei 2017.10.20
};

//解析Eclipse格式文件暂定接口 by 孙亮 2017.11.16 from RIPED
struct Eclipse_Data_Interface
{
	//节点
	vector<float> X;
	vector<float> Y;
	vector<float> Z;
	//属性
	map<string, vector<float> > Pros;
};
Eclipse_Data_Interface Read_Eclipse(string fn);

//************************************
// 函数名:  OnReadShiWenData
// 功能:    读取网格模型数据
// 参数:    模型文件名strFileName,读取到值存放到modeldata中 
// 返回值:  成功 true,失败 false
// 说明:    
// 日期:    2017-08-16
// 作者:    刘少兵
//************************************
//DECLARE_GMEDLL_EXPORT_IMPORT bool OnReadShiWenData(const string& strFileName, iPeraModData &modeldata);
bool OnReadShiWenData(const string& strFileName, iPeraModData &modeldata);

//************************************
// 函数名:  OnWriteShiWenData
// 功能:    保存网格模型数据
// 参数:    模型文件名strFileName,modeldata中的数据存储到磁盘文件
// 返回值:  成功 true,失败 false
// 说明:      
// 日期:    2017-08-16
// 作者:    刘少兵
//************************************
//DECLARE_GMEDLL_EXPORT_IMPORT bool OnWriteShiWenData(const string& strFileName, iPeraModData &modeldata);
bool OnWriteShiWenData(const string& strFileName, iPeraModData &modeldata);