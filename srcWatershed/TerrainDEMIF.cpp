#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <fstream>
using std::ifstream;

#include<zlib.h>
#include "TerrainDEMIF.H"


#include "NamespaceHeader.H"

TerrainDEMIF::TerrainDEMIF(const ProblemDomain&     a_domain,
			   const int                a_interpType,
			   const RealVect&          a_finestDx,
			   const RealVect&          a_coarseDx,
			   const std::string&       a_demFile,
			   const Real               a_bottomElevation,
			   const RealVect           a_origin,
			   const Real               a_truncElev,
			   const Real               a_highGround,
			   const Real               a_verticalScale,
			   const Real               a_cellvalue,
			   const int                a_edgeType)
{
    CH_TIME("TerrainDEMIF::TerrainDEMIF");
    m_highGround    = a_highGround;
    m_verticalScale = a_verticalScale;
    m_origin =  a_origin;
    m_cellvalue = a_cellvalue;
    
    m_ncell         = a_domain.size();
    m_finestDx      = a_finestDx;
    m_coarseDx      = a_coarseDx;
    m_interpType    = a_interpType;
    m_DEM=NULL;
    m_edgeType      = a_edgeType;

    if (m_interpType==1)
    {
	m_doCubic = false;
    }
    else
    {
	m_doCubic = true;
    }
    
    m_bottomElevation  = a_bottomElevation;
    m_truncElev     = a_truncElev;
    

    // Read in the whole file
    readDEM(a_demFile);

    //fill in stray "wet" cells, change nodata value to land
    // and convert DEM to meters.
    fixDEM();
}

TerrainDEMIF::TerrainDEMIF(const TerrainDEMIF& a_inputIF)
{
    CH_TIME("TerrainDEMIF::TerrainDEMIF");
  m_ncols         = a_inputIF.m_ncols;
  m_nrows         = a_inputIF.m_nrows;
  m_xllcorner     = a_inputIF.m_xllcorner;
  m_yllcorner     = a_inputIF.m_yllcorner;
  m_hx            = a_inputIF.m_hx;
  m_hy            = a_inputIF.m_hy;
  m_cellvalue     = a_inputIF.m_cellvalue;
  m_finestDx      = a_inputIF.m_finestDx;
  m_coarseDx      = a_inputIF.m_coarseDx;
  m_ncell         = a_inputIF.m_ncell;
  m_interpType    = a_inputIF.m_interpType;
  m_highGround    = a_inputIF.m_highGround;
  m_minDEM        = a_inputIF.m_minDEM;
  m_maxDEM        = a_inputIF.m_maxDEM;
  m_bottomElevation = a_inputIF.m_bottomElevation;
  m_truncElev     = a_inputIF.m_truncElev;
  m_verticalScale = a_inputIF.m_verticalScale;
  m_edgeType      = a_inputIF.m_edgeType;

  if (m_interpType==1)
  {
      m_doCubic = false;
  }
  else
  {
      m_doCubic = true;
  }

  //allocate and copy memory for DEM matrix
  m_DEM = new Real *[m_ncols];

  for (int i = 0; i < m_ncols;i++)
  {
      m_DEM[i] = new Real[m_nrows];
      for (int j = 0; j < m_nrows;j++)
      {
          m_DEM[i][j] = a_inputIF.m_DEM[i][j];
      }
  }
  m_hx2    = pow(m_hx,2);
  m_hx3    = pow(m_hx,3);
  m_hy2    = pow(m_hy,2);
  m_hy3    = pow(m_hy,3);
}

TerrainDEMIF::~TerrainDEMIF()
{
    CH_TIME("TerrainDEMIF::~TerrainDEMIF");

    for (int i = 0; i < m_ncols;i++)
    {
	delete [] m_DEM[i];
    }
    delete [] m_DEM;
}

//Reads in a Digital Elevation Model file
bool TerrainDEMIF::readDEM(const std::string& a_demFile)
{
    CH_TIME("TerrainDEMIF::readDEM");
    bool fileread;
    //pout()<<"Reading DEM Header...\n";
    ifstream DEM_file;
    gzFile zip_file_stream=NULL;

    //check if mobility file is arc type
    bool zip_flag=false;
    std::string file_ext = a_demFile.substr(a_demFile.find_last_of(".") + 1);
    if(file_ext == "gz") 
    {
	zip_flag=true;
	zip_file_stream=gzopen(a_demFile.c_str(), "r");
	if (!zip_file_stream) MayDay::Error("TerrainDEMIF:: Could not open zipped file.");
	gzread(zip_file_stream, &m_ncols, sizeof(int));
	gzread(zip_file_stream, &m_nrows, sizeof(int));
	gzread(zip_file_stream, &m_hx, sizeof(Real));
	gzread(zip_file_stream, &m_hy, sizeof(Real));
    }
    else
    {
	DEM_file.open(a_demFile.c_str()); //the input file
	if (!DEM_file.good())
	{
	    MayDay::Abort("Bad DEM_file in readDEM");
	}
	char astring[1024];
	DEM_file>>astring;
	DEM_file>>m_ncols;
	DEM_file>>m_nrows;
	DEM_file>>astring;
	DEM_file>>m_hx;
	DEM_file>>m_hy;
    }

    // pout()<<m_hx<<"\t"<<m_hy<<endl;
    // pout()<<m_ncols<<"\t"<<m_nrows<<endl; 

//allocate memory for DEM matrix
    m_DEM = new Real *[m_ncols];
    for (int i = 0; i < m_ncols;i++)
    {
	m_DEM[i] = new Real[m_nrows];
    }
    
    m_hx2    = pow(m_hx,2);
    m_hx3    = pow(m_hx,3);
    m_hy2    = pow(m_hy,2);
    m_hy3    = pow(m_hy,3);
    
    //shrink the DEM domain, so we will
    //use additional points for the interpolation
    int sz=0;
    // if (m_interpType==1)
    // {
    // 	sz=0;
    // }
    // else
    // {
    // 	sz=0;
    // }
    m_xllcorner = m_origin[0]-(Real)sz*m_hx;
    m_yllcorner = m_origin[1]-(Real)sz*m_hy;
    
    
    
  // if (SpaceDim==2)
  // {
  //     if (m_nrows!=1)
  //     {//2D column data...switch m_nrows and m_ncols
  //         int temp = m_nrows;
  //         m_nrows = m_ncols;
  //         m_ncols = temp;
  //     }
  //     else if (m_ncols==1)
  //     {
  //         MayDay::Error("TerrainDEMIF::ncols and nrows are 1...need a bigger dataset");
  //     }
  //     CH_assert(m_nrows==1);
  // }
  // else if (SpaceDim==3 && m_nrows==1)
  // {
  //     MayDay::Error("TerrainDEMIF::nrows must not equal 1 for 3D");
  // }
  
  //make sure the grid is big enough
  RealVect domainHi = m_ncell;
  domainHi *=m_finestDx;
  domainHi +=m_origin;
  
  if ((domainHi[0] > (m_xllcorner+(Real)(m_ncols-sz)*m_hx)) ||
      (domainHi[1] > (m_yllcorner+(Real)(m_nrows-sz)*m_hy)))
  {
      MayDay::Error("TerrainDEMIF::DEM file does not provide a sufficient grid size for interpolation.");
  }

  {
      m_minDEM =  1e20;
      m_maxDEM = -1e20;
      //load in the data
      if (zip_flag)
      { 
	  unsigned int bytes_read=0;
	  for (int i = 0; i < m_ncols;i++)
	  {
	      bytes_read +=gzread(zip_file_stream, m_DEM[i], sizeof(double)*m_nrows);
	  }

	  if (bytes_read!=sizeof(double)*m_nrows*m_ncols)  MayDay::Error("TerrainDEMIF::Error reading gz file.");
	  
	  for (int j=0;j<m_nrows;j++)
	  {
	      for (int i=0;i<m_ncols;i++)
	      {
		  
		  //vertical shift and scale
		  Real val = (-m_bottomElevation + m_DEM[i][j]*m_cellvalue)*m_verticalScale;
		  m_DEM[i][j] = val;
		  if (val < m_minDEM)
		  {
		      m_minDEM = val; //this is the highest spot so far...
		  }
		  if (val > m_maxDEM)
		  {
		      m_maxDEM = val; //this is the deepest spot so far...
		  }
	      }
	  }
      }
      else
      {
	  for (int j=0;j<m_nrows;j++)
	  {
	      for (int i=0;i<m_ncols;i++)
	      {
		  
		  DEM_file>>m_DEM[i][j];
		  //pout()<<i<<","<<j<<"\t"<<m_DEM[i][j];
		  //vertical shift and scale
		  Real val = (-m_bottomElevation + m_DEM[i][j]*m_cellvalue)*m_verticalScale;
		  // pout()<<"\t"<<val<<endl;
		  
		  m_DEM[i][j] = val;
		  if (val < m_minDEM)
		  {
		      m_minDEM = val; //this is the highest spot so far...
		  }
		  if (val > m_maxDEM)
		  {
		      m_maxDEM = val; //this is the deepest spot so far...
		  }
	      }
	  }
      }
      pout()<<"  DEM min = "<<m_minDEM<<endl;
      pout()<<"  DEM max = "<<m_maxDEM<<endl;
      
  }
  //load in the data
  if (zip_flag)
  {
      gzclose(zip_file_stream);
  }
  else
  {
      DEM_file.close();
  }

  fileread = true;
  return fileread;
}

// change nodata value to land
// and convert DEM to meters.
bool TerrainDEMIF::fixDEM()
{
    CH_TIME("TerrainDEMIF::fixDEM");
  bool DEMfixed = true;

  int kt = 0;
  for (int j=0;j<m_nrows;j++)
  {
      for (int i=0;i<m_ncols;i++)
      {
          if (m_DEM[i][j] > m_truncElev)
	  {
              m_DEM[i][j] = m_highGround;
              kt++;
	  }
      }
  }
  pout()<<"Filled: "<<kt<<" shallow cells."<<endl;
  return DEMfixed;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

Real TerrainDEMIF::value(const RealVect& a_point) const
{
    CH_TIME("TerrainDEMIF::value(point)");

  Real retval =0.;

  RealVect x = a_point;
  //we need to convert it from reference frame with unit grid spacing to physical frame
  x *=m_coarseDx;

  //where xllcorner is the x-lower-left corner of the cell-centered DEM
  //find nearest cell in DEM (recall that the DEM values are cell centered)
  Real xi,yj;
  xi = (x[0]-m_xllcorner)/m_hx - 0.5;
  yj = (x[1]-m_yllcorner)/m_hy - 0.5;
  Real xiFloor = floor(xi);
  Real yjFloor = floor(yj);

  int i,j;
  // centered on i+1/2
  i = (int)xiFloor;
  j = (int)yjFloor;

  int sz;//this is the size of the interp stencil
  if (m_doCubic)
  {//cubic
      sz = 1;
  }
  else
  {//linear
      sz = 0;
  }

  //the following checks to make sure we are inside of the DEM database 
  if (i+sz<0 || i+sz>m_ncols || j+sz<0 || j+sz>m_nrows) 
  {
       MayDay::Error("TerrainDEMIF::DEM file does not provide a sufficient grid size for interpolation_.");
  }
  else
  {
      if (m_interpType==1)
      { //bilinear interpolation
          RealVect a   = RealVect::Zero;
          Real xLo = m_xllcorner + (i+0.5)*m_hx; //this is the x position of the low side of this dem cell (i,j)
          Real yLo = m_yllcorner + (j+0.5)*m_hy; //this is the y position of the low side of this dem cell (i,j)
	  
          //a is the (bi)linear coefficient
          a[0] = (x[0]-xLo)/m_hx;
          a[1] = (x[1]-yLo)/m_hy;

	  const Real fLoLo = getValueFromDEM(i,  j);
	  const Real fHiLo = getValueFromDEM(i+1,j);
	  
	  Real rLo = (1.0-a[0])*fLoLo + a[0]*fHiLo;
	  retval = rLo;
	  
          const Real fLoHi = getValueFromDEM(i,j+1);
          const Real fHiHi = getValueFromDEM(i+1,j+1);
          Real rHi = (1.0-a[0])*fLoHi + a[0]*fHiHi;
          retval   = (1.0-a[1])*rLo   + a[1]*rHi;
      }
      else if (m_interpType==2)
      { //bicubic interpolation type1 (this performs cubic interpolation in 1D, matching all values)
          //now let's do some bicubic interpolation...

	  const Real f1 = getValueFromDEM(i-1,j-1);
          const Real f2 = getValueFromDEM(i,  j-1);
          const Real f3 = getValueFromDEM(i+1,j-1);
          const Real f4 = getValueFromDEM(i+2,j-1);

          const Real f5 = getValueFromDEM(i-1,j);
          const Real f6 = getValueFromDEM(i,  j);
          const Real f7 = getValueFromDEM(i+1,j);
          const Real f8 = getValueFromDEM(i+2,j);
          const Real f9 = getValueFromDEM(i-1,j+1);
          const Real f10= getValueFromDEM(i,  j+1);
          const Real f11= getValueFromDEM(i+1,j+1);
          const Real f12= getValueFromDEM(i+2,j+1);

          const Real f13= getValueFromDEM(i-1,j+2);
          const Real f14= getValueFromDEM(i,  j+2);
          const Real f15= getValueFromDEM(i+1,j+2);
          const Real f16= getValueFromDEM(i+2,j+2);

	  Real xd=x[0]-m_xllcorner-m_hx*(float(i+1)); //this is the x-distance from (i+1/2,j+1/2)
          Real x2 = pow(xd,2);
          Real x3 = pow(xd,3);

          Real yd=x[1]-m_yllcorner-m_hy*(float(j+1)); //this is the y-distance from (i+1/2,j+1/2)
          Real y2 = pow(yd,2);
          Real y3 = pow(yd,3);

          //cubic type1 in x direction: for row j-1
          Real r1 = (-3.*m_hx3*( f1 - 9.*(f2 + f3) + f4)    +
                     2.*m_hx2*( f1 -27.*(f2 - f3) - f4)*xd  +
                     12.*m_hx*( f1 -    (f2 + f3) + f4)*x2 +
                     8.*(-f1 + 3.*(f2 - f3) + f4)*x3)/(48.*m_hx3);

          //cubic type1 in x direction: for row j
          Real r2 = (-3.*m_hx3*( f5 - 9.*(f6 + f7) + f8)    +
                     2.*m_hx2*( f5 -27.*(f6 - f7) - f8)*xd  +
                     12.*m_hx*( f5 -    (f6 + f7) + f8)*x2 +
                     8.*(-f5 + 3.*(f6 - f7) + f8)*x3)/(48.*m_hx3);
          retval = r2;

          //cubic type1 in x direction: for row j+1
          Real r3 = (-3.*m_hx3*( f9 - 9.*(f10 + f11) + f12)    +
                     2.*m_hx2*( f9 -27.*(f10 - f11) - f12)*xd  +
                     12.*m_hx*( f9 -    (f10 + f11) + f12)*x2 +
                     8.*(-f9 + 3.*(f10 - f11) + f12)*x3)/(48.*m_hx3);

          //cubic type1 in x direction: for row j+2
          Real r4 = (-3.*m_hx3*( f13 - 9.*(f14 + f15) + f16)    +
                     2.*m_hx2*( f13 -27.*(f14 - f15) - f16)*xd  +
                     12.*m_hx*( f13 -    (f14 + f15) + f16)*x2 +
                     8.*(-f13 + 3.*(f14 - f15) + f16)*x3)/(48.*m_hx3);

          //cubic type1 in y direction!!!
          retval  = (-3.*m_hy3*( r1 - 9.*(r2 + r3) + r4)    +
                     2.*m_hy2*( r1 -27.*(r2 - r3) - r4)*yd  +
                     12.*m_hy*( r1 -    (r2 + r3) + r4)*y2 +
                     8.*(-r1 + 3.*(r2 - r3) + r4)*y3)/(48.*m_hy3);
      }
      else if (m_interpType==3)
      { //bicubic interpolation type2 (this performs cubic interpolation in 1D, matching values and derivatives)

	  const Real f1 = getValueFromDEM(i-1,j-1);
          const Real f2 = getValueFromDEM(i,  j-1);
          const Real f3 = getValueFromDEM(i+1,j-1);
          const Real f4 = getValueFromDEM(i+2,j-1);
          const Real f5 = getValueFromDEM(i-1,j);
          const Real f6 = getValueFromDEM(i,  j);
          const Real f7 = getValueFromDEM(i+1,j);
          const Real f8 = getValueFromDEM(i+2,j);
          const Real f9 = getValueFromDEM(i-1,j+1);
          const Real f10= getValueFromDEM(i,  j+1);
          const Real f11= getValueFromDEM(i+1,j+1);
          const Real f12= getValueFromDEM(i+2,j+1);
          const Real f13= getValueFromDEM(i-1,j+2);
          const Real f14= getValueFromDEM(i,  j+2);
          const Real f15= getValueFromDEM(i+1,j+2);
          const Real f16= getValueFromDEM(i+2,j+2);


          Real xd=x[0]-m_xllcorner-m_hx*(float(i+1)); //this is the x-distance from (i+1/2,j+1/2)
          Real x2 = pow(xd,2);
          Real x3 = pow(xd,3);
          Real yd=x[1]-m_yllcorner-m_hy*(float(j+1)); //this is the y-distance from (i+1/2,j+1/2)
          Real y2 = pow(yd,2);
          Real y3 = pow(yd,3);

          //cubic type2 in x direction: for row j-1
          Real r1 = (-f1 + 9.*(f2 + f3) - f4)   /    16.  +
            ( f1 -11.*(f2 - f3) - f4)*xd /(m_hx *8.) +
            ( f1 -    (f2 + f3) + f4)*x2/(m_hx2*4.) +
            (-f1 + 3.*(f2 - f3) + f4)*x3/(m_hx3*2.);
          //cubic type2 in x direction: for row j
          Real r2 = (-f5 + 9.*(f6 + f7) - f8)   /    16.  +
            ( f5 -11.*(f6 - f7) - f8)*xd /(m_hx *8.) +
            ( f5 -    (f6 + f7) + f8)*x2/(m_hx2*4.) +
            (-f5 + 3.*(f6 - f7) + f8)*x3/(m_hx3*2.);
          retval = r2;
          //cubic type2 in x direction: for row j+1
          Real r3 = (-f9 + 9.*(f10 + f11) - f12)   /    16.  +
            ( f9 -11.*(f10 - f11) - f12)*xd /(m_hx *8.) +
            ( f9 -    (f10 + f11) + f12)*x2/(m_hx2*4.) +
            (-f9 + 3.*(f10 - f11) + f12)*x3/(m_hx3*2.);

          //cubic type2 in x direction: for row j+2
          Real r4 = (-f13 + 9.*(f14 + f15) - f16)   /    16.  +
            ( f13 -11.*(f14 - f15) - f16)*xd /(m_hx *8.) +
            ( f13 -    (f14 + f15) + f16)*x2/(m_hx2*4.) +
            (-f13 + 3.*(f14 - f15) + f16)*x3/(m_hx3*2.);

          //cubic type2 in y direction!!!
          retval  = (-r1 + 9.*(r2 + r3) - r4)   /    16.  +
            ( r1 -11.*(r2 - r3) - r4)*yd /(m_hy *8.) +
            ( r1 -    (r2 + r3) + r4)*y2/(m_hy2*4.) +
            (-r1 + 3.*(r2 - r3) + r4)*y3/(m_hy3*2.);
      }
      else
      {
          MayDay::Error("Bad InterpType in TerrainDEMIF");
      }
  }

  //pout()<<a_point<<"\t retval="<<retval;
  retval /=m_coarseDx[2];
  retval -= a_point[2];//*m_coarseDx[2];
  
  //pout()<<"\t retval_2="<<-retval<<endl;
  
  return -retval;
}

Real
TerrainDEMIF::value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
		    const IndexTM<Real,GLOBALDIM>& a_point) const
{
    CH_TIME("TerrainDEMIF::value");
  Real retval= LARGEREALVAL;
  int derivativeOrder = a_partialDerivative.sum();
  
  if (derivativeOrder == 0)
  {
      retval = value(a_point);
  }
  else if (derivativeOrder > 1)
  {
      //derivative function is not implemented
      //however we need this function to be defined in IFData class
      //setting it to zero
      
      //pout()<<"derivative function is not implemented yet"<<endl;
      retval = 0.0;
  }
  
  // Change the sign to change inside to outside
  // if (!m_inside && derivativeOrder > 0)
  // {
  //   retval = -retval;
  // }

  return retval;
}

BaseIF* TerrainDEMIF::newImplicitFunction() const
{
  TerrainDEMIF* demPtr = new TerrainDEMIF(*this);

  return static_cast<BaseIF*>(demPtr);
}


Real TerrainDEMIF::getValueFromDEM(const int a_i, const int a_j) const
{
    int i,j;
    if (a_i<0)
    {
	if(m_edgeType==1)
	{
	    i=abs(a_i);
	}
	else
	{
	    i=0;
	}
	
    }else if (a_i>m_ncols-1)
    {
	if(m_edgeType==1)
	{
	    i=2*m_ncols-a_i-2;
	}
	else
	{
	    i=m_ncols-1;
	}
	
    }else
    {
	i=a_i;
    }

    if (a_j<0)
    {
	if(m_edgeType==1)
	{
	    j=abs(a_j);
	}
	else
	{
	    j=0;
	}
	
    }else if (a_j>m_nrows-1)
    {
	if(m_edgeType==1)
	{
	    j=2*m_nrows-a_j-2;
	}
	else
	{
	    j=m_nrows-1;
	}
	
    }else
	
    {
	j=a_j;
    }
    

    return m_DEM[i][j];
}


#include "NamespaceFooter.H"
