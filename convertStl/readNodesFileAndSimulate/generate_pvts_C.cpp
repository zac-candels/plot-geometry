#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <stdio.h>
#include <sys/stat.h>
using namespace std;

int outtype,numfluid,numsolute,numtemp;
 
const char* dataSp = "./OutData_PV";
const char* dataSp1= "/PV_";
const char* dataSpF1 = "/F_PV_";
const char* dataSp2= "/PCM_"; 
const char* dataF = "./File"; 
const char* dataF2= "/File_Info";
const char* dataF3= "/File_Ext"; 
 
const char* dataSp3= "./PVData";
 
const int kki = 0, TMax = 1.0e8;
const int Step1 = 1e4, Step2 = 5e4, Step3 = 5e4;
const int MaxS1 = 1e5, MaxS2 = 1e6;
const int ns=1e4, nEs=1e6; // initial iteration
const int nEs1=5e2, nEs2=1e3; // middle iteration
 
class cGatherData
{
    // creates parallel VTK files (lbtout%.4d.pvts) to link files created by
    // individual processes (lbout*.vts), where %.4d is the saved time step
    // the data structure follows the VTK XML standard for structured grids
  
    public:

    cGatherData()
    {
        getsize();
        gatherallVTK();
    }

    ~cGatherData()
    {

    }

    int getsize();
    int gatherallVTK();
  
protected:
 
    int sizeofsys;
    int sizeofver;
    int ntx, nty, ntz;
    int *xs, *xe, *ys, *ye, *zs, *ze;
};

int cGatherData::getsize()
{
    // this procedure is required to avoid mismatches caused by different machines 
    char buf[80];
    char issue[22]; 
    int x1, x2, y1, y2, z1, z2;
    int j=0;
    int n=0;
    int Num=0;

    if (outtype<0 || outtype>4) 
    {
        cout<<"Which output type was produced?"<<endl;
        cout<<"(0 = all, 1 = density, 2 = solid charge, 3 = liquid charge, 4 = Vos, 5 = solid)"<<endl;
        while (outtype<0 || outtype>4) cin >> outtype; 
    }
  
    //  determine size of pieces 
    sprintf(buf, "%s%s", dataF, dataF2);  
    ifstream myfile(buf); 
    if(!myfile)
    {
        cout<<"buf "<<buf<<" \n"; 
        cout<<"error opening "<<"File_info"<<" file\n"; 
        exit(1);
    } 
    while (!myfile.eof() && j<2)
    {
        myfile >> issue >> x1;

        if(!strcmp(issue, "sizeofSystem")) 
        {
            sizeofsys = int(x1);
            j++;
        } 
    }
    myfile.close();
 
    /////////////////
    xs=new int[sizeofsys];
    xe=new int[sizeofsys];
    ys=new int[sizeofsys];
    ye=new int[sizeofsys];
    zs=new int[sizeofsys];
    ze=new int[sizeofsys];
     
    //  determine extent of each piece and entire system 
    sprintf(buf, "%s%s", dataF, dataF3);  
    ifstream myfile1(buf);  
    if(!myfile1)
    {
      cout<<"error opening "<<"File_ext"<<" file\n"; 
      exit(1);
    }
    ntx = 0; nty = 0; ntz = 0;
    while (!myfile1.eof())
    {
        myfile1 >> issue >> x1 >> x2 >> y1 >> y2 >> z1 >> z2;
        
        for(int ip=0; ip<sizeofsys; ip=ip+1) 
        {
            sprintf(buf, "extent_%d", ip);
            if(!strcmp(issue, buf)) 
            {
                xs[ip] = x1;  xe[ip] = x2;
                ys[ip] = y1;  ye[ip] = y2;
                zs[ip] = z1;  ze[ip] = z2;
                if(x2>ntx) ntx = x2;
                if(y2>nty) nty = y2;
                if(z2>ntz) ntz = z2;
            } 
        }
    }
    myfile1.close();
 
    //  determine number of VTK-files
    n=0; Num=0;
    while(1)
    {  
        sizeofver = Num;
        
        sprintf(buf, "%s%s%d%s%d", dataSp, dataSp1, n, dataSp2, 0);
        ifstream myfile(buf); 
        if(!myfile) break;
        else
        {
            Num=Num+1;  
            n=n+Step1;
        } 
    }
 
    // deal with unequal numbers of grid points among processes
    struct stat results;
    n=0; Num=0;
    for(Num=0; Num<sizeofver; Num=Num+1)
    { 
        for(int ip=0; ip<sizeofsys; ip=ip+1)
        {
            sprintf(buf, "%s%s%d%s%d", dataSp, dataSp1, n, dataSp2, ip);
 
            if(stat(buf, &results) != 0) 
            {
                cout<<"error opening "<<buf<<" file \n";
                exit(1);
            }
        }
        n=n+Step1;
    }
 
    return 0;
}

int cGatherData::gatherallVTK()
{ 
    // write .pvts file for each timestep
    char namebuf[80];
 
    ostringstream command;
    command<<"mkdir -p "<<dataSp3;
    system(command.str().c_str());
 
    int n=0;
    for(int Num=0; Num<sizeofver; Num=Num+1)
    {
        sprintf(namebuf, "%s%s%.8d.pvts", dataSp3, dataSp2, n);
        ofstream ofile(namebuf);

        ofile<<"<?xml version=\"1.0\"?>"<<endl;
        //ofile<<"<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">"<<endl;
        //ofile<<"<PStructuredGrid WholeExtent=\"0 "<<ntx<<" 0 "<<nty<<" 0 "<<ntz<<"\" GhostLevel=\"1\">"<<endl;
        ofile<<"<VTKFile type=\"PRectilinearGrid\" version=\"0.1\" byte_order=\"BigEndian\">"<<endl;
        ofile<<"<PRectilinearGrid WholeExtent=\"0 "<<ntx<<" 0 "<<nty<<" 0 "<<ntz<<"\" GhostLevel=\"1\">"<<endl;
        ofile<<"<PPointData Scalars=\"phase_field\" Vectors=\"velocity\">"<<endl;

        if (outtype==0) 
        {
            // for(int iprop=0; iprop<numfluid; iprop++)
            // ofile<<"<PDataArray Name=\"density\" type=\"Float32\"/>"<<endl;
            // for(int iprop=0; iprop<numfluid; iprop++)
            // ofile<<"<PDataArray Name=\"fraction_"<<iprop<<"\" type=\"Float32\"/>"<<endl;
            // for(int iprop=0; iprop<numsolute; iprop++)
            // ofile<<"<PDataArray Name=\"concentration_"<<iprop<<"\" type=\"Float32\"/>"<<endl;
            // if(numtemp==1)
            // ofile<<"<PDataArray Name=\"temperature\" type=\"Float32\"/>"<<endl;

            ofile<<"<PDataArray Name=\"VOF\" type=\"Float32\"/>"<<endl; 
            ofile<<"<PDataArray Name=\"Solid1\" type=\"Float32\"/>"<<endl;
            ofile<<"<PDataArray Name=\"Solid2\" type=\"Float32\"/>"<<endl;
            ofile<<"<PDataArray Name=\"u\" type=\"Float32\" NumberOfComponents=\"3\"/>"<<endl;
        }
        else if (outtype==1) ofile<<"<PDataArray Name=\"density\" type=\"Float32\"/>"<<endl;
        else if (outtype==2) ofile<<"<PDataArray Name=\"solid charge\" type=\"Float32\"/>"<<endl;
        else if (outtype==3) ofile<<"<PDataArray Name=\"liquid charge\" type=\"Float32\"/>"<<endl;
        else if (outtype==4) ofile<<"<PDataArray Name=\"Vos\" type=\"Float32\"/>"<<endl;
        else if (outtype==5) ofile<<"<PDataArray Name=\"solid\" type=\"Float32\"/>"<<endl;
        else ofile<<"<PDataArray Name=\"temperature\" type=\"Float32\"/>"<<endl;

        // ofile<<"<PDataArray Name=\"velocity\" type=\"Float32\" NumberOfComponents=\"3\"/>"<<endl;
        // ofile<<"<PDataArray Name=\"phase_field\" type=\"Int32\"/>"<<endl;
        ofile<<"</PPointData>"<<endl;
        // ofile<<"<PPoints>"<<endl;
        // ofile<<"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>"<<endl;
        // ofile<<"</PPoints>"<<endl;

        ofile<<"<PCoordinates>"<<endl;
        ofile<<"<PDataArray type=\"Float32\"/>"<<endl;
        ofile<<"<PDataArray type=\"Float32\"/>"<<endl;
        ofile<<"<PDataArray type=\"Float32\"/>"<<endl;
        ofile<<"</PCoordinates>"<<endl;
        
        for(int i=0; i<sizeofsys; i=i+1) 
        {
            ofile<<"<Piece Extent=\""<<xs[i]<<" "<<xe[i]<<" "<<ys[i]<<" "<<ye[i]<<" "<<zs[i]<<" "<<ze[i]
                <<"\" Source=\"."<<dataSp<<dataSp1<<n<<"/"<<dataSp2<<i<<"\"/>"<<endl;
        }

        //    ofile<<"</PStructuredGrid>"<<endl;

        ofile<<"</PRectilinearGrid>"<<endl;
        ofile<<"</VTKFile>"<<endl;
        ofile.close();

        //////////////
        n=n+Step1;
    }

  return 0;
}

int main(int argc, char* argv[])
{
    outtype=0;

    if (argc>1) outtype = atoi(argv[1]);

    cGatherData aa;

    return 0;
}
