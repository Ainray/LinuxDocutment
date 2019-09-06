////////////////////////////////////////////////////////////////
// 版权所有 (c) 2017, 石文软件有限公司
// 保留所有权利
// 
// 文件名称：ShiWen_iPeraMod_DataInterface.cpp
// 摘    要：
//			盆模系统网格模型接口实现
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
// 最后修改：2017-09-2f6
////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "ShiWen_iPeraMod_DataInterface.h"

//**************第一部分：基本功能函数******************//
//功能：申请内存
template <class T>
static inline bool GetBuffer(T * &pDataBuffer, size_t nSize)
{
	if (pDataBuffer)
	{
		delete[]pDataBuffer;
		pDataBuffer = NULL;
	}
	try
	{
		pDataBuffer = new T[nSize];
		return true;
	}
	catch (...)
	{
		return false;
	}
}

//功能：释放内存
template <class T>
static inline void ReleaseBuffer(T * &pDataBuffer)
{
	if (pDataBuffer != NULL)
	{
		delete[]pDataBuffer;
	}
	pDataBuffer = NULL;
}

//功能：获取属性名称
string GetProName(char *proName, int nSize)
{
	int pos = -1;
	for (int i = 0; i < nSize; i++)
	{
		char ch = proName[i];
		if (ch == '\0')
		{
			pos = i;
			break;
		}
	}
	string strName;
	if (pos != -1)
	{
		strName.assign(proName, pos);
	}
	return strName;
}

//功能：保存属性名称
static void SaveProName(const string& strName, ofstream &file, int nProNameLength)
{
	int nLength = (int)strName.length();
	char *name = new char[nProNameLength];
	memcpy(name, strName.data(), nLength);
	name[nLength] = '\0';
	file.write(name, nProNameLength * sizeof(char));
	delete[]name;
	name = NULL;
}
//***********************************************************************//

//***********第二部分：从文件读石文数据到iPeraMod接口*******************//
//步骤：01
//功能：读取文件

//步骤：02 （优化方案 by 孙亮）
//功能：读取石文数据接口数据头
//输入：dataFile文件
//返回：石文数据接口
ShiWenFileHead ReadShiWenFileHead(ifstream& dataFile)
{
	ShiWenFileHead outHead;
	dataFile.read((char*)&outHead, sizeof(ShiWenFileHead));
	return outHead;
}

//步骤：03
//功能：读取属性名称
//输入：dataFile文件,nNameLength名称固定长度,nProNum属性个数
//返回：所有属性名称
std::vector<string> ReadShiWenFileProNames(ifstream& dataFile, int nNameLength, int nProNum)
{
	char *proNames = new char[nNameLength * nProNum];
	dataFile.read(proNames, nNameLength * nProNum * sizeof(char));
	std::vector< string > names(nProNum);
	for (int i = 0; i < nProNum; i++)
	{
		int offset = i*nNameLength;
		names[i] = GetProName(proNames + offset, nNameLength);
	}
	delete[]proNames;
	proNames = NULL;
	return names;
}

//步骤：04
//功能：读取每一大层的小层层数
//输入：dataFile文件,LayerNum个数
//返回：每层小层数
std::vector<int> ReadShiWenFileLayerNumbers(ifstream& dataFile, int LayerNum)
{
	std::vector<int> vecEachLayerNumber(LayerNum);
	dataFile.read((char*)&vecEachLayerNumber[0], LayerNum * sizeof(int));
	return vecEachLayerNumber;
}

//功能：读取网格模型各层的名称
//输入：dataFile文件,nNameLength名称固定长度,nLayerNum层数
//返回：网格模型的层名称 
//add by xuwei 2017.10.20
std::vector<std::string> ReadShiWenFileLayerNames(std::ifstream& dataFile, int nNameLength, int nLayerNum)
{
	char *layerNames = new char[nNameLength * nLayerNum];
	dataFile.read(layerNames, nNameLength * nLayerNum * sizeof(char));
	std::vector<std::string> names(nLayerNum);
	for (int i = 0; i < nLayerNum; i++)
	{
		int offset = i*nNameLength;
		names[i] = GetProName(layerNames + offset, nNameLength);
	}
	delete[]layerNames;
	layerNames = NULL;
	return names;
}

//步骤：05 (优化方案 by 孙亮)
//功能：读取属性数据
//输入：dataFile文件, nPros 属性个数，nNodes节点个数
//输出: 属性数组
//vector<vector<float> > ReadShiWenFileProData(ifstream& dataFile, size_t nPros, int nNodes)
vector<float> ReadShiWenFileProData(ifstream& dataFile,size_t nPros, int nNodes)
{
#ifdef _TWODIMENSION_
	vector<vector<float> > Pros;
	vector<float> data(nNodes);	
	for (size_t i = 0; i < nPros; i++)
		{
			dataFile.read((char*)&data[0], nNodes * sizeof(float));
			Pros.push_back(data);
		}
#else
	vector<float> Pros(nNodes*nPros);
	dataFile.read((char *)&Pros[0],nNodes*nPros*sizeof(float));
#endif
	return Pros;
}

//步骤：06
//功能：关闭文件


//************************************
// 步骤：00 （优化方案）
// 函数名:  OnReadShiWenData 
// 功能:    读取石文网格模型数据
// 参数:    模型文件名strFileName,读取到值存放到modeldata中 
// 返回值:  成功 true,失败 false
// 说明:    
// 日期:    2017-09-22
// 作者:    孙亮
//************************************
//DECLARE_GMEDLL_EXPORT_IMPORT bool OnReadShiWenData(const std::string& strFileName, iPeraModData &modeldata)
bool OnReadShiWenData(const std::string& strFileName, iPeraModData &modeldata)
{
	//第一步：打开文件
	ifstream dataFile(strFileName, ios::binary);
	if (!dataFile.is_open())
	{
		cout << "opening input file failed !" << endl;
		return false;
	}
	//第二步：读取文件头
	ShiWenFileHead outHead = ReadShiWenFileHead(dataFile);
	if (outHead.nFileFlag != 0x8001)
	{
		dataFile.close();
		cout << "input file formet is wrong !" << endl;
		return false;
	}
	//第三步:读取属性名称
	vector<string> proNames = ReadShiWenFileProNames(dataFile, outHead.nProNameLength, outHead.nProNum);
	//第四步:读取层信息
	vector<int> vecEachLayerNumber = ReadShiWenFileLayerNumbers(dataFile, outHead.nLayerNum);
	std::vector<std::string> vecLayerNames = ReadShiWenFileLayerNames(dataFile, outHead.nProNameLength, outHead.nLayerNum);//读取层名 by xuwei 2017.10.20
	//第五步：读取属性信息
	int nNodes = outHead.nXNodeNum * outHead.nYNodeNum * outHead.nZNodeNum;
//	vector<vector<float> > proValues = ReadShiWenFileProData(dataFile, outHead.nProNum, nNodes);
	vector<float> proValues = ReadShiWenFileProData(dataFile, outHead.nProNum, nNodes);

	//第六步：关闭文件
	dataFile.close();
	//第七步：传递返回值
	modeldata.ShiWenFileHead = outHead;
	modeldata.proNames = proNames;
	modeldata.proValues = proValues;
	modeldata.vecEachLayerNumber = vecEachLayerNumber;
	modeldata.vecEachLayerName = vecLayerNames;
	return true;
}
//***************************************************************//

//***********第三部分：写石文数据到文件*******************//

//步骤：01
//功能：打开文件


//步骤：02 (优化方案)
//功能：写模型数据头
//by 孙亮

bool WriteShiWenFileHead(ofstream& dataFile, ShiWenFileHead &outHead)
{
	try
	{
		dataFile.write((char*)&outHead, sizeof(ShiWenFileHead));
		return true;
	}
	catch (...)
	{
		return false;
	}
}

//步骤：03 （优化方案）
//功能：写属性名称
//输入：dataFile文件,proNames字符数组,nNameLength名称固定长度
//by 孙亮
bool WriteShiWenFileProNames(ofstream& dataFile, vector<string> proNames, int nProNameLength)
{
	try
	{
		for (unsigned int i = 0; i<proNames.size(); i++)
		{
			SaveProName(proNames[i], dataFile, nProNameLength);
		}
		return true;
	}
	catch (...)
	{
		return false;
	}
}

//步骤：04 (优化方案)
//功能：写每一大层的小层层数
//输入：dataFile文件，vecEachLayerNumber每层layer个数
//by 孙亮
bool WriteShiWenFileLayerNumbers(ofstream& dataFile, vector<int> vecEachLayerNumber)
{
	try
	{
		bool bRet = false;
		if (vecEachLayerNumber.size() > 0)
		{
			dataFile.write((char*)&vecEachLayerNumber[0], vecEachLayerNumber.size() * sizeof(int));
			bRet = true;
		}
		return bRet;
	}
	catch (...)
	{
		return false;
	}
}

//功能：写层名称
//输入：dataFile文件,layerNames字符数组,nNameLength名称固定长度
//add by xuwei 2017.10.20
bool WriteShiWenFileLayerNames(std::ofstream& dataFile, std::vector<std::string>& layerNames, int nProNameLength)
{
	try
	{
		for (unsigned int i = 0; i < layerNames.size(); i++)
		{
			SaveProName(layerNames[i], dataFile, nProNameLength);
		}
		return true;
	}
	catch (...)
	{
		return false;
	}
}


//步骤：05 （优化方案）
//功能：写属性数据
//输入：dataFile文件,vector<vector<float> > data属性数组
//输出:成功true ,失败false
//by 孙亮
#ifdef _TWODIMENSION_
bool WriteShiWenFileProData(ofstream& dataFile, vector<vector<float> > proValues)
{
	try
	{
	 	for (size_t i = 0; i < proValues.size(); i++)
		{
			std::vector<float> data = proValues[i];
			dataFile.write((char*)&data[0], data.size() * sizeof(float));
		}
		return true;
	}
	catch (...)
	{
		return false;
	}
}
#else
bool WriteShiWenFileProData(ofstream& dataFile, vector<float> proValues)
{
	try
	{
		dataFile.write((char*)&proValues[0], proValues.size() * sizeof(float));
		return true;
	}
	catch (...)
	{
		return false;
	}
}
#endif
//步骤：06
//功能：关闭文件

//************************************
// 步骤：00(优化方案)
// 函数名:  OnWriteShiWenData
// 功能:    保存网格模型数据
// 参数:    模型文件名strFileName,modeldata中的数据存储到磁盘文件
// 返回值:  成功 true,失败 false
// 说明:    
// 日期:    2017-8-16
// 作者:    刘少兵, 孙亮20170922
//************************************
//DECLARE_GMEDLL_EXPORT_IMPORT bool OnWriteShiWenData(const std::string& strFileName, iPeraModData &modeldata)
bool OnWriteShiWenData(const std::string& strFileName, iPeraModData &modeldata)
{
	//第一步：打开文件
	ofstream dataFile(strFileName, ios::binary);
	if (!dataFile.is_open()) 
	{
		cout << "opening output file failed !" << endl;
		return false;
	}

	//参数检验
	if (modeldata.vecEachLayerNumber.size() == 0
		|| modeldata.ShiWenFileHead.nProNum != (int)modeldata.proNames.size())
//		|| modeldata.ShiWenFileHead.nProNum != (int)modeldata.proValues.size())
	{
		cout << "parameter checking failed !" << endl;
		return false;
	}

	//第二步：写文件头
	if (!WriteShiWenFileHead(dataFile, modeldata.ShiWenFileHead))
	{
		cout << "writing file head failed !" << endl;
		return false;
	}
		
	//第三步: 写属性名称
	if (!WriteShiWenFileProNames(dataFile, modeldata.proNames, modeldata.ShiWenFileHead.nProNameLength))
	{
		cout << "writing pro names failed !" << endl;
		return false;
	}

	//第四步: 写各大层层数信息
	if (!WriteShiWenFileLayerNumbers(dataFile, modeldata.vecEachLayerNumber))
	{
		cout << "writing each layer numbers failed !" << endl;
		return false;
	}

	//~~~: 写各大层层名, by xuwei 2017.10.20
	if (!WriteShiWenFileLayerNames(dataFile, modeldata.vecEachLayerName, modeldata.ShiWenFileHead.nProNameLength))
	{
		cout << "writing each layer names failed !" << endl;
		return false;
	}

	//第五步：写属性数据
	if(!WriteShiWenFileProData(dataFile, modeldata.proValues))
	{
		cout << "writing pro values failed !" << endl;
		return false;
	}

	//第六步：关闭文件
	dataFile.close();
	cout << "Shiwen Grid Data is written successfully !" << endl;
	return true;
}
//********************************************************************//

//解析Eclipse格式文件暂定接口代码 by 孙亮 2017.11.16 from RIPED
Eclipse_Data_Interface Read_Eclipse(string fn)
{
	Eclipse_Data_Interface Model;

	int _mI, _mJ, _mK;
	std::string kw, tmp;

	std::ifstream fin(fn);
	if (!fin) {
		printf("Error opening file %s", fn.c_str());
		//return false;
	}

	//-------------------read dimensions-----------------------//
	while (1) {
		std::getline(fin, kw);
		if (0 == kw.compare(0, 8, "SPECGRID"))
			break;
	}
	fin >> _mI >> _mJ >> _mK;
	printf("grid dimension is: %d  %d  %d\n", _mI, _mJ, _mK);

	float *pDepth = (float*)malloc(sizeof(float)*_mI * 2);
	float *pAttrib = (float*)malloc(sizeof(float)*_mI); ;
	//---------------------------------------------------------//


	//-------------------parse ZCOORN section-------------------//
	printf("begin parse ZCOORN section\n");
	printf("--------------------------\n");
	while (1) {
		std::getline(fin, kw);
		if (0 == kw.compare(0, 5, "ZCORN"))
			break;
	}

	float *p;
	float *a = (float*)malloc(sizeof(float)*_mI * 2);
	p = &a[0];
	for (int i = 0; i < _mI; i++) {
		*p++ =(float) (10000 * i);   // x interval -- 100000
		*p++ =(float)( 10000 * (1 + i));
	}
	float *b = (float*)malloc(sizeof(float)*_mJ * 2);
	p = &b[0];
	for (int i = 0; i < _mJ; i++) {
		*p++ =(float) (10000 * i);   // y interval -- 100000
		*p++ =(float)( 10000 * (1 + i));
	}

	int NodeID = 0;
	//vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New(); //----- vtk point setting --------------//
	for (int k = 0; k < _mK * 2; ++k) {
		for (int j = 0; j < _mJ * 2; ++j) {
			float *p = pDepth;
			for (int m = 0; m<2 * _mI / 6; m++) {
				fin >> p[6 * m] >> p[6 * m + 1] >> p[6 * m + 2] >> p[6 * m + 3] >> p[6 * m + 4] >> p[6 * m + 5];
			}
			std::getline(fin, tmp);

			for (int i = 0; i < _mI * 2; ++i) {
				float pt[] = { a[i], b[j], p[i] };
				//points->InsertNextPoint(pt); //----- vtk point setting --------------//

				//----- model node setting --------------// 
				{
					Model.X.push_back(a[i]);
					Model.X.push_back(b[i]);
					Model.X.push_back(p[i]);
				}
				//----- model node setting --------------// 
			}
		}
	}

	// Add the points and hexahedron to an unstructured grid
	//vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	//uGrid->SetPoints(points);
	//---------------------------------------------------------//

	//--------------setting cells -------------------//
	int _cid = 0;  //cell id
				   //point id of each hexa
	for (int k = 0; k<_mK; k++) {
		for (int j = 0; j<_mJ; j++) {
			for (int i = 0; i<_mI; i++) {
				// Create a hexahedron from the points
				//vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();

				/*
				*        unsigned int idx00 = 22080*(2*k)+120*(2*j)+(2*i);
				*        unsigned int idx01 = 22080*(2*k)+120*(2*j)+(2*i+1);
				*        unsigned int idx02 = 22080*(2*k)+120*(2*j+1)+(2*i+1);
				*        unsigned int idx03 = 22080*(2*k)+120*(2*j+1)+(2*i);
				*
				*        unsigned int idx10 = 22080*(2*k+1)+120*(2*j)+(2*i);
				*        unsigned int idx11 = 22080*(2*k+1)+120*(2*j)+(2*i+1);
				*        unsigned int idx12 = 22080*(2*k+1)+120*(2*j+1)+(2*i+1);
				*        unsigned int idx13 = 22080*(2*k+1)+120*(2*j+1)+(2*i);
				*
				*/

				int PN = _mI*_mJ * 4;

				unsigned int idx00 = PN*(2 * k) + 2 * _mI*(2 * j) + (2 * i);
				unsigned int idx01 = PN*(2 * k) + 2 * _mI*(2 * j) + (2 * i + 1);
				unsigned int idx02 = PN*(2 * k) + 2 * _mI*(2 * j + 1) + (2 * i + 1);
				unsigned int idx03 = PN*(2 * k) + 2 * _mI*(2 * j + 1) + (2 * i);

				unsigned int idx10 = PN*(2 * k + 1) + 2 * _mI*(2 * j) + (2 * i);
				unsigned int idx11 = PN*(2 * k + 1) + 2 * _mI*(2 * j) + (2 * i + 1);
				unsigned int idx12 = PN*(2 * k + 1) + 2 * _mI*(2 * j + 1) + (2 * i + 1);
				unsigned int idx13 = PN*(2 * k + 1) + 2 * _mI*(2 * j + 1) + (2 * i);

				//----- vtk cell setting --------------//
				//hex->GetPointIds()->SetId(0, idx00); hex->GetPointIds()->SetId(1, idx01);
				//hex->GetPointIds()->SetId(2, idx02); hex->GetPointIds()->SetId(3, idx03);

				//hex->GetPointIds()->SetId(4, idx10); hex->GetPointIds()->SetId(5, idx11);
				//hex->GetPointIds()->SetId(6, idx12); hex->GetPointIds()->SetId(7, idx13);

				//uGrid->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
				//----- vtk cell setting --------------//

				//----- model cell setting --------------// 
				//{
					//Cell cell;
					//cell.id = _cid; cell.Type = "hexahedron";

					//cell.Ids.push_back(idx00); cell.Ids.push_back(idx01);
					//cell.Ids.push_back(idx02); cell.Ids.push_back(idx03);

					//cell.Ids.push_back(idx10); cell.Ids.push_back(idx11);
					//cell.Ids.push_back(idx12); cell.Ids.push_back(idx13);

					//Cells.push_back(cell);
					//_cid++;
				//}
				//----- model node setting --------------// 

			}
		}
	}
	//---------------------------------------------------------//

	//-----------------read props-----------------------------//
	printf("begin parse GRAIN1 section\n");
	{
		std::string tmp;
		while (1) {
			std::getline(fin, kw);
			if (0 == kw.compare(0, 6, "GRAIN1"))
				break;
		}
		float *p = pAttrib;
		//vtkSmartPointer<vtkFloatArray> prop = vtkSmartPointer<vtkFloatArray>::New();
		//prop->SetName("GRAIN1");
		//vtkCellData* celldat = uGrid->GetCellData();

		for (int k = 0; k < _mK; ++k) {
			for (int j = 0; j < _mJ; ++j) {
				for (int m = 0; m<_mI / 6; m++)
					fin >> p[6 * m] >> p[6 * m + 1] >> p[6 * m + 2] >> p[6 * m + 3] >> p[6 * m + 4] >> p[6 * m + 5];
				std::getline(fin, tmp);

				for (int i = 0; i < _mI; ++i)
				{
					if (p[i]<-9000)
						p[i] = NAN;

					//prop->InsertNextTuple1(p[i]); //----- vtk prop setting --------------//

					//----- grid attrib value
					{
						Model.Pros["GRAIN1"].push_back(p[i]);
					}
					//----- grid attrib value
				}
			}
		}
		//celldat->AddArray(prop);
	}

	printf("begin parse GRAIN2 section\n");
	{
		//=====================parse GRAIN2 section
		std::string tmp;
		while (1) {
			std::getline(fin, kw);
			if (0 == kw.compare(0, 6, "GRAIN2"))
				break;
		}
		float *p = pAttrib;
		//vtkSmartPointer<vtkFloatArray> prop = vtkSmartPointer<vtkFloatArray>::New();
		//prop->SetName("GRAIN2");
		//vtkCellData* celldat = uGrid->GetCellData();

		for (int k = 0; k < _mK; ++k) {
			for (int j = 0; j < _mJ; ++j) {
				for (int m = 0; m<_mI / 6; m++)
					fin >> p[6 * m] >> p[6 * m + 1] >> p[6 * m + 2] >> p[6 * m + 3] >> p[6 * m + 4] >> p[6 * m + 5];
				std::getline(fin, tmp);

				for (int i = 0; i < _mI; ++i)
				{
					if (p[i]<-9000)
						p[i] = NAN;

					//prop->InsertNextTuple1(p[i]); //----- vtk prop setting --------------//

					//----- grid attrib value
					{
						Model.Pros["GRAIN2"].push_back(p[i]);
					}
					//----- grid attrib value
				}
			}
		}
		//celldat->AddArray(prop);
	}

	printf("begin parse PORO section\n");
	{
		std::string tmp;
		while (1) {
			std::getline(fin, kw);
			if (0 == kw.compare(0, 4, "PORO"))
				break;
		}

		float *p = pAttrib;

		//vtkSmartPointer<vtkFloatArray> poro = vtkSmartPointer<vtkFloatArray>::New();
		//poro->SetName("PORO");
		//vtkCellData* celldat = uGrid->GetCellData();

		for (int k = 0; k < _mK; ++k) {
			for (int j = 0; j < _mJ; ++j) {
				for (int m = 0; m<_mI / 6; m++)
					fin >> p[6 * m] >> p[6 * m + 1] >> p[6 * m + 2] >> p[6 * m + 3] >> p[6 * m + 4] >> p[6 * m + 5];
				std::getline(fin, tmp);

				for (int i = 0; i < _mI; ++i)
				{
					if (p[i]<-9000)
						p[i] = NAN;

					//poro->InsertNextTuple1(p[i]); //----- vtk prop setting --------------//

					//----- grid attrib value
					{
						Model.Pros["PORO"].push_back(p[i]);
					}
					//----- grid attrib value
				}
			}
		}
		//celldat->AddArray(poro);
	}

	printf("begin parse PREM section\n");
	{
		std::string tmp;
		while (1) {
			std::getline(fin, kw);
			if (0 == kw.compare(0, 4, "PERM"))
				break;
		}
		float *p = pAttrib;
		//vtkSmartPointer<vtkFloatArray> perm = vtkSmartPointer<vtkFloatArray>::New();
		//perm->SetName("PERM");
		//vtkCellData* celldat = uGrid->GetCellData();

		for (int k = 0; k < _mK; ++k) {
			for (int j = 0; j < _mJ; ++j) {
				for (int m = 0; m<_mI / 6; m++)
					fin >> p[6 * m] >> p[6 * m + 1] >> p[6 * m + 2] >> p[6 * m + 3] >> p[6 * m + 4] >> p[6 * m + 5];
				std::getline(fin, tmp);

				for (int i = 0; i < _mI; ++i)
				{
					if (p[i]<-9000)
						p[i] = NAN;

					//perm->InsertNextTuple1(p[i]); //----- vtk prop setting --------------//

					//----- grid attrib value
					{
						Model.Pros["PERM"].push_back(p[i]);
					}
					//----- grid attrib value
				}
			}
		}
		//celldat->AddArray(perm);
	}

	// ---- grid statics
	//nNodes = Nodes.size();
	//nCells = Cells.size();
	//nNPros = NPros.size();

	// Write the file
	//vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	//writer->SetFileName("strata.vtu");
	//writer->SetInputData(uGrid);
	//writer->Write();

	return Model;
}