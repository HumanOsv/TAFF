#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
//#include <vector>

#include <sys/types.h>
#include <sys/stat.h>


#include "FileManager.h"
using namespace std;

class approach {

	string BossOne[5];	//struct with charge info from orbitals. //change Koop, orbitalweighted and fd, 2 if per fucntion, fd is new
	struct charges{				
		string charge;
		double n_charge;
		charges *next, *prev;
	};
	charges *root, *homo_position, *lumo_position;
	float Che_pot, TotalEnergy;	//chemical potential and  Total energy.
	int count;
	public:
		
		approach();//constructor de la clase...
		void Koopmans(string file, string name, int,int);
		void OrbitalWeighted(string , float, string, int);		
		void Finite_Differences(string, string, string,string,int);;
		void GetEnergy(ifstream &file, int);
		void GetHOMO_LUMO(ifstream &file, int);
		void GetListCharges(ifstream &file, int, string);
		void Show_HOMO_LUMO_Info();
		void ChemicalPotential();
		void GlobalHardness();
		void ElectrophilicityIndex(float );
		void Electronegativity();
		void GlobalSoftness(float );
		float IonizationPotencial();
		float ElectroAffinity();
		void GAP();
		float ElectrophilicityDonator(float, float);
		float ElectrophilicityAceptor(float, float);
		void ElectrophilicityNet(float, float);
		float IonizationPotencialFD(float , float );
		float ElectroAffinityFD(float , float );
		void ChemicalPotentialFD(float , float );
		void GlobalHardnessFD(float , float );
		float getTotalEnergy();
		void Directory(string);
		void Density(string, string, string, string, string, fileManager*);
		void Get_Dual(string,string,string);
		
};//end class definition.
approach::approach(){	//obligatory constructor.
	root=new charges;
}
void approach::Directory(string route){	//create a folder if it doesnt exist.
	struct stat st = {0};

	if (stat(route.c_str(), &st) == -1) {
    	mkdir(route.c_str(), 0777);	

	}
}
void approach::Density(string name, string file, string cube, string app, string type, fileManager *printer){
	string order="5\n1\n8\n"+cube+".cube\n2\nx\n";
	string name_rhoCube="f"+type+app+"-"+name+"-rho.cube";
	printer->Multiwfn_Pipeline(order,file, "tmp.txt");
	rename("density.cub",name_rhoCube.c_str());
//estoy cambiando la orden del ulti abajo, parte 2, falta cambiar el -3 por el rho.cube
	order="17\n1\n2\nx\n";	
	if(printer->Multiwfn_Pipeline(order,cube+".cube", "tmp_d"+type+".txt")==1){
		printer->DeleteFile("tmp_d"+type+".txt");
		printer->Multiwfn_Pipeline("17\n1\n2\n1\n2\n-1\n"+name_rhoCube+"\nx\n",cube+".cube", "tmp_d"+type+".txt");
	}else{
		printer->DeleteFile("tmp_d"+type+".txt");		
		printer->Multiwfn_Pipeline("17\n1\n2\n2\n-1\n"+name_rhoCube+"\nx\n",cube+".cube", "tmp_d"+type+".txt");	
	}
}

void approach::Koopmans(string file, string name, int work_option, int dual){
	//vector<fileManager> fm;

	//fileManager *printer=new fileManager();	//.cube names.
	fileManager *printer=new fileManager();
	//printer->MultiwfnRoute();
	string order, fminus="f-koop-"+name+".cube", fplus="f+koop-"+name+".cube", fzero="f0rad_koop-"+name+".cube";
	//.cibe pipeline f-
	order="6\n26\n0\n0.0\n"+BossOne[0]+"\n1.0\n00\n-1\n5\n1\n3\n2\nx\n";	
	printer->Multiwfn_Pipeline(order,file, "tmp.txt");
	rename("density.cub",fminus.c_str());
	//f+
	order="6\n26\n0\n0.0\n"+BossOne[1]+"\n1.0\n00\n-1\n5\n1\n3\n2\nx\n";
	printer->Multiwfn_Pipeline(order,file, "tmp.txt");
	rename("density.cub",fplus.c_str());
	//
	if(dual==1){
		Get_Dual(fplus,fminus,name);
	}
	//f0rad
	order="13\n11\n12\n"+fplus+"\n0\n"+fzero+"\nx\n";
	printer->Multiwfn_Pipeline(order,fminus, "tmp.txt");
	//get the basin integra, volumen and coordinates
	if(work_option==1 || work_option==3){
		order="17\n1\n2\nx\n";	
		if(printer->Multiwfn_Pipeline(order,fminus, "tmp_f-.txt")==1){
			printer->DeleteFile("tmp_f-.txt");
			printer->Multiwfn_Pipeline("17\n1\n2\n1\n2\n0\n-3\nx\n",fminus, "tmp_f-.txt");
		}else{
			printer->DeleteFile("tmp_f-.txt");		
			printer->Multiwfn_Pipeline("17\n1\n2\n2\n0\n-3\nx\n",fminus, "tmp_f-.txt");	
		}
		if(printer->Multiwfn_Pipeline(order,fplus, "tmp_f+.txt")==1){
			printer->DeleteFile("tmp_f+.txt");
			printer->Multiwfn_Pipeline("17\n1\n2\n1\n2\n0\n-3\nx\n",fplus, "tmp_f+.txt");
		}else{
			printer->DeleteFile("tmp_f+.txt");
			printer->Multiwfn_Pipeline("17\n1\n2\n2\n0\n-3\nx\n",fplus, "tmp_f+.txt");
		}
		if(printer->Multiwfn_Pipeline(order,fzero, "tmp_f0.txt")==1){
			printer->DeleteFile("tmp_f0.txt");
			printer->Multiwfn_Pipeline("17\n1\n2\n1\n2\n0\n-3\nx\n",fzero, "tmp_f0.txt");
		}else{
			printer->DeleteFile("tmp_f0.txt");
			printer->Multiwfn_Pipeline("17\n1\n2\n2\n0\n-3\nx\n",fzero, "tmp_f0.txt");	
		}
		//integral txt data
		printer->Extrac_N_Write_Integral_DAT("tmp_f-.txt", "f-koop-"+name,"Koopmans_cube");
		printer->Extrac_N_Write_Integral_DAT("tmp_f+.txt", "f+koop-"+name,"Koopmans_cube");
		printer->Extrac_N_Write_Integral_DAT("tmp_f0.txt", "f0rad_koop-"+name,"Koopmans_cube");
	}
	//wrtite -integral.txt and vmd file
	//here route to cube

	fminus="f-koop-"+name;
	fplus="f+koop-"+name;
	fzero="f0rad_koop-"+name;
	if(work_option==2 || work_option==3){
		//Density order with Multi
		Density(name, file,fminus,"koop","-",printer);
		Density(name, file,fplus,"koop","+",printer);
		Density(name, file,fzero,"koop","0rad_",printer);
		printer->Density_Data_File("tmp_d-.txt",fminus,"Koopmans_cube");
		printer->Density_Data_File("tmp_d+.txt",fplus,"Koopmans_cube");
		printer->Density_Data_File("tmp_d0rad_.txt",fzero,"Koopmans_cube");
		
	}
	//detele tmp files.
	printer->DeleteFile("tmp_f+.txt");
	printer->DeleteFile("tmp_f-.txt");
	printer->DeleteFile("tmp_f0.txt");
	printer->DeleteFile("tmp_d+.txt");
	printer->DeleteFile("tmp_d-.txt");
	printer->DeleteFile("tmp_d0rad_.txt");
	printer->DeleteFile("tmp.txt");
	
	//////LINUX//////////////
	Directory("Koopmans_cube");
	Directory("Koopmans_txt");
	//move files to folders.
	rename(string(fminus+".cube").c_str(),string("Koopmans_cube/f-koop-"+name+".cube").c_str());
	rename(string(fplus+".cube").c_str(),string("Koopmans_cube/f+koop-"+name+".cube").c_str());
	rename(string(fzero+".cube").c_str(),string("Koopmans_cube/f0rad_koop-"+name+".cube").c_str());

	rename(string(fminus+"-VMD.vmd").c_str(),string("Koopmans_cube/f-koop-"+name+"-VMD.vmd").c_str());
	rename(string(fplus+"-VMD.vmd").c_str(),string("Koopmans_cube/f+koop-"+name+"-VMD.vmd").c_str());
	rename(string(fzero+"-VMD.vmd").c_str(),string("Koopmans_cube/f0rad_koop-"+name+"-VMD.vmd").c_str());

	rename(string(fminus+"-integral.txt").c_str(),string("Koopmans_txt/f-koop-"+name+"-integral.txt").c_str());
	rename(string(fplus+"-integral.txt").c_str(),string("Koopmans_txt/f+koop-"+name+"-integral.txt").c_str());
	rename(string(fzero+"-integral.txt").c_str(),string("Koopmans_txt/f0rad_koop-"+name+"-integral.txt").c_str());

	//Density folders
	rename(string(fminus+"-rho.cube").c_str(),string("Koopmans_cube/f-koop-"+name+"-rho.cube").c_str());
	rename(string(fplus+"-rho.cube").c_str(),string("Koopmans_cube/f+koop-"+name+"-rho.cube").c_str());
	rename(string(fzero+"-rho.cube").c_str(),string("Koopmans_cube/f0rad_koop-"+name+"-rho.cube").c_str());

	rename(string(fminus+"-rho-VMD.vmd").c_str(),string("Koopmans_cube/f-koop-"+name+"-rho-VMD.vmd").c_str());
	rename(string(fplus+"-rho-VMD.vmd").c_str(),string("Koopmans_cube/f+koop-"+name+"-rho-VMD.vmd").c_str());
	rename(string(fzero+"-rho-VMD.vmd").c_str(),string("Koopmans_cube/f0rad_koop-"+name+"-rho-VMD.vmd").c_str());

	rename(string(fminus+"-rho.txt").c_str(),string("Koopmans_txt/f-koop-"+name+"-rho.txt").c_str());
	rename(string(fplus+"-rho.txt").c_str(),string("Koopmans_txt/f+koop-"+name+"-rho.txt").c_str());
	rename(string(fzero+"-rho.txt").c_str(),string("Koopmans_txt/f0rad_koop-"+name+"-rho.txt").c_str());


}
void approach::OrbitalWeighted(string file, float gaussian, string name, int work_option){
	int HOMO_position_charge= atoi(BossOne[0].c_str()), aux=0;	//variables used
	charges *helper=root;
	long double sum_fMinus=0, sum_fPlus=0,sus,division,square,expo;	//variables from equations.
	long double OccupiedOWValues[HOMO_position_charge], UnoccupiedOWValues[count -HOMO_position_charge];
	fileManager *printer=new fileManager();
	ostringstream gaus;	//to add gaussian value to the name of the .cube, dat, vmd files.
	gaus << gaussian;
	name=name+"_"+gaus.str()+"_";

//	cout << "total: "<<count<<" HOMO "<< HOMO_position_charge<< " Else "<< count - HOMO_position_charge<<endl;
	string order1, order2, fminus="f-orbitalweighted-"+name+".cube", fplus="f+orbitalweighted-"+name+".cube", fzero="f0rad_orbitalweighted-"+name+".cube";
	order1="6\n26\n";
	order2=order1;
	if ( helper != 0 ) { //Makes sure there is a place to start
	  	while ( helper->next != 0 ) {
	  		ostringstream convert, convertEXP;
	  		aux++;							//orbitalweighted equation exponencial(-((chemicalpotencial-charge)**2)/gaussian))
	  		if(aux<=HOMO_position_charge){	//sum equation until HOMO.
		  		sus=Che_pot-helper->n_charge;
		  		square=sus*sus;
		  		division=square/gaussian;
		  		expo= expl(division*(-1.0));
	  			OccupiedOWValues[aux-1]=expo;
	  			convert << aux ;
	  			convertEXP << expo;
	  			sum_fMinus=sum_fMinus+expo;	//sumatory
	  			order1=order1+convert.str()+"\n"+convertEXP.str()+"\n";	//add to Multiwfn pipeline
	  			order2=order2+convert.str()+"\n0.0\n";
	  		}else{	//sum equation LUMO and forwards
	  			//cout <<"LUMO"<<endl;
		  		sus=Che_pot-helper->n_charge;
		  		square=sus*sus;
		  		division=square/gaussian;
	  
	  			expo= expl(division*(-1.0));
	  			if(division<=40){
	  				UnoccupiedOWValues[aux-HOMO_position_charge-1]= expo;
	  				//order2=order2+convert.str()+"\n"+convertEXP.str()+"\n";
	  				sum_fPlus=sum_fPlus+expo;
	  			}else{
	  				break;
	  			}
	  		}
	  		helper = helper->next;
	  	}
	}
	ostringstream  convertEXP;
	convertEXP << sum_fMinus;
	order1=order1+"00\n-1\n5\n1\n3\n6\n"+convertEXP.str()+"\n2\nx\n";//Multiwfn pipeline for real orbitals.
	///
	for (int i = 0; i <= (aux-HOMO_position_charge-1); ++i)
	{
		stringstream aux, value;
		aux << (i+HOMO_position_charge+1);
		value << (UnoccupiedOWValues[i]/sum_fPlus);
		order2=order2+aux.str()+"\n"+value.str()+"\n";
	}
	///
	//Here we make the array with the contriibutions of each orbital.
	//
	
	order2=order2+"00\n-1\n5\n1\n3\n2\nx\n"; 
	//pipeline to .cube files f-
	printer->Multiwfn_Pipeline(order1,file, "tmp.txt");
	rename("density.cub",fminus.c_str());
	//f+
	//printer->Multiwfn_Pipeline(order2,file, "tmp.txt");
	//rename("density.cub",fplus.c_str());
	
	if(work_option==1 || work_option==3){
		order1=	"17\n1\n2\nx\n";
		//order1="17\n1\n2\n2\n0\n-3\nx\n";
		// printer->Multiwfn_Pipeline(order1,fminus, "tmp_f-.txt"); //basin, integral, volumen and coordinates values.	
		if(printer->Multiwfn_Pipeline(order1,fminus, "tmp_f-.txt")==1){
			printer->DeleteFile("tmp_f-.txt");
			printer->Multiwfn_Pipeline("17\n1\n2\n1\n2\n0\n-3\nx\n",fminus, "tmp_f-.txt");
		}else{
			printer->DeleteFile("tmp_f-.txt");
			printer->Multiwfn_Pipeline("17\n1\n2\n2\n0\n-3\nx\n",fminus, "tmp_f-.txt");
		}/*
		if(printer->Multiwfn_Pipeline(order1,fplus, "tmp_f+.txt")==1){
			printer->DeleteFile("tmp_f+.txt");
			printer->Multiwfn_Pipeline("17\n1\n2\n1\n2\n0\n-3\nx\n",fplus, "tmp_f+.txt");
		}else{
			printer->DeleteFile("tmp_f+.txt");
			printer->Multiwfn_Pipeline("17\n1\n2\n2\n0\n-3\nx\n",fplus, "tmp_f+.txt");
		}*/
		//write -integral.txt and vmd files.
		printer->Extrac_N_Write_Integral_DAT("tmp_f-.txt", "f-orbitalweighted-"+name,"OrbitalWeighted_cube");	//files creations
		//printer->Extrac_N_Write_Integral_DAT("tmp_f+.txt", "f+orbitalweighted-"+name,"OrbitalWeighted_cube");	//files creations

		printer->PrintPcenOW(OccupiedOWValues, HOMO_position_charge,UnoccupiedOWValues, (aux-HOMO_position_charge-1), sum_fMinus, sum_fPlus,name );
		//printer->Extrac_N_Write_Integral_DAT("tmp_f0.txt", "f0rad_orbitalweighted-"+name,"OrbitalWeighted_cube");	//files creations
	}
	fminus="f-orbitalweighted-"+name;
//	fplus="f+orbitalweighted-"+name;
	if(work_option==2 || work_option==3){
		Density(name, file,fminus,"orbitalweighted","-",printer);
		printer->Density_Data_File("tmp_d-.txt",fminus,"OrbitalWeighted_txt_cube");
	}
	
	//delete tmp files.
	printer->DeleteFile("tmp_f+.txt");
	printer->DeleteFile("tmp_f-.txt");
	printer->DeleteFile("tmp_f0.txt");
	printer->DeleteFile("tmp.txt");
	printer->DeleteFile("tmp_d-.txt");
	//LINUX
	Directory("OrbitalWeighted_cube");
	Directory("OrbitalWeighted_txt");
	//move files to folders
	rename(string(fminus+".cube").c_str(),string("OrbitalWeighted_cube/f-orbitalweighted-"+name+".cube").c_str());
	//rename(string(fplus+".cube").c_str(),string("OrbitalWeighted_cube/f+orbitalweighted-"+name+".cube").c_str());
	//rename(string(fzero+".cube").c_str(),string("OrbitalWeighted_cube/f0rad_orbitalweighted-"+name+".cube").c_str());

	rename(string(fminus+"-VMD.vmd").c_str(),string("OrbitalWeighted_cube/f-orbitalweighted-"+name+"-VMD.vmd").c_str());
	//rename(string(fplus+"-VMD.vmd").c_str(),string("OrbitalWeighted_cube/f+orbitalweighted-"+name+"-VMD.vmd").c_str());
	//rename(string(fzero+"-VMD.vmd").c_str(),string("OrbitalWeighted_cube/f0rad_orbitalweighted-"+name+"-VMD.vmd").c_str());

	rename(string(fminus+"-integral.txt").c_str(),string("OrbitalWeighted_txt/f-orbitalweighted-"+name+"-integral.txt").c_str());
	//rename(string(fplus+"-integral.txt").c_str(),string("OrbitalWeighted_txt/f+orbitalweighted-"+name+"-integral.txt").c_str());
	//rename(string(fzero+"-integral.txt").c_str(),string("OrbitalWeighted_txt/f0rad_orbitalweighted-"+name+"-integral.txt").c_str());

	rename(string(fminus+"-rho.cube").c_str(),string("OrbitalWeighted_cube/f-orbitalweighted-"+name+"-rho.cube").c_str());

	//rename(string(fminus+"-rho-VMD.vmd").c_str(),string("OrbitalWeighted_cube/f-orbitalweighted-"+name+"-rho-VMD.vmd").c_str());

	//rename(string(fminus+"-rho.txt").c_str(),string("OrbitalWeighted_txt/f-orbitalweighted-"+name+"-rho.txt").c_str());

}
void approach::Finite_Differences(string fileZero, string fileMinus, string filePlus, string Realname, int work_option){
	//vector<fileManager> fm;
	string order, File_name_MINUScube="Anion.cube", File_name_PLUScube="Cation.cube", File_name_ZEROcube=Realname+".cube";
	//fileManager *printer=new fileManager();
	fileManager *printer=new fileManager();
	//printer->MultiwfnRoute();
	//neutral.cube pipeline
	order="5\n1\n3\n2\nx\n";
	printer->Multiwfn_Pipeline(order,fileZero, "tmp.txt");
	rename("density.cub",File_name_ZEROcube.c_str());
	if(filePlus!="0."){
		//cation.cube pipeline
		order="5\n1\n8\n"+File_name_ZEROcube+"\n2\nx\n";
		printer->Multiwfn_Pipeline(order,filePlus, "tmp.txt");
		rename("density.cub",File_name_PLUScube.c_str());
		//f- pipeline
		order="13\n11\n4\n"+File_name_PLUScube+"\n0\nf-_dif-"+File_name_ZEROcube+"\nx\n";
		printer->Multiwfn_Pipeline(order,File_name_ZEROcube, "tmp.txt");
	}if(fileMinus!="0."){
		//anion.cube pipeline
		order="5\n1\n8\n"+File_name_ZEROcube+"\n2\nx\n";
		printer->Multiwfn_Pipeline(order,fileMinus, "tmp.txt");
		rename("density.cub",File_name_MINUScube.c_str());
		//f+ pipeline.
		order="13\n11\n4\n"+File_name_ZEROcube+"\n0\nf+_dif-"+File_name_ZEROcube+"\nx\n";
		printer->Multiwfn_Pipeline(order,File_name_MINUScube, "tmp.txt");
	}
	printer->DeleteFile(File_name_MINUScube);	//delete neutral cation anion cube
	printer->DeleteFile(File_name_ZEROcube);
	printer->DeleteFile(File_name_PLUScube);
	if(fileMinus!="0." && filePlus!="0."){	//forad
		order="13\n11\n12\nf+_dif-"+File_name_ZEROcube+"\n0\nf0rad_dif-"+File_name_ZEROcube+"\nx\n";
		printer->Multiwfn_Pipeline(order,"f-_dif-"+File_name_ZEROcube, "tmp.txt");	
	}
	if(work_option==1 || work_option==3){
	//basin, integral, volumen and coordinates values.
	order="17\n1\n2\nx\n";
	//	order="17\n1\n2\n2\n0\n-3\nx\n";
		if(filePlus!="0."){
			if(printer->Multiwfn_Pipeline(order,"f-_dif-"+File_name_ZEROcube, "tmp_f-.txt")==1){
				printer->DeleteFile("tmp_f-.txt");
				printer->Multiwfn_Pipeline("17\n1\n2\n1\n2\n0\n-3\nx\n","f-_dif-"+File_name_ZEROcube, "tmp_f-.txt");
			}else{
				printer->DeleteFile("tmp_f-.txt");
				printer->Multiwfn_Pipeline("17\n1\n2\n2\n0\n-3\nx\n","f-_dif-"+File_name_ZEROcube, "tmp_f-.txt");					
			}
			printer->Extrac_N_Write_Integral_DAT("tmp_f-.txt", "f-_dif-"+Realname,"FiniteDifferences_cube");
		}if(fileMinus!="0."){
			if(printer->Multiwfn_Pipeline(order,"f+_dif-"+File_name_ZEROcube, "tmp_f+.txt")==1){
				printer->DeleteFile("tmp_f+.txt");
				printer->Multiwfn_Pipeline("17\n1\n2\n1\n2\n0\n-3\nx\n","f+_dif-"+File_name_ZEROcube, "tmp_f+.txt");			
			}else{
				printer->DeleteFile("tmp_f+.txt");
				printer->Multiwfn_Pipeline("17\n1\n2\n2\n0\n-3\nx\n","f+_dif-"+File_name_ZEROcube, "tmp_f+.txt");					
			}
			printer->Extrac_N_Write_Integral_DAT("tmp_f+.txt", "f+_dif-"+Realname,"FiniteDifferences_cube");
		}if(fileMinus!="0." && filePlus!="0."){
			if(printer->Multiwfn_Pipeline(order,"f0rad_dif-"+File_name_ZEROcube, "tmp_f0.txt")==1){
				printer->DeleteFile("tmp_f0.txt");
				printer->Multiwfn_Pipeline("17\n1\n2\n1\n2\n0\n-3\nx\n","f0rad_dif-"+File_name_ZEROcube, "tmp_f0.txt");
			}else{
				printer->DeleteFile("tmp_f0.txt");
				printer->Multiwfn_Pipeline("17\n1\n2\n2\n0\n-3\nx\n","f0rad_dif-"+File_name_ZEROcube, "tmp_f0.txt");					
			}
			printer->Extrac_N_Write_Integral_DAT("tmp_f0.txt", "f0rad_dif-"+Realname,"FiniteDifferences_cube");
		}
	}
	File_name_MINUScube="f-_dif-"+Realname;
	File_name_PLUScube="f+_dif-"+Realname;
	File_name_ZEROcube="f0rad_dif-"+Realname;
	if(work_option==2 || work_option==3){
		if(filePlus!="0."){
			Density(Realname, fileMinus,File_name_MINUScube,"_dif","-",printer);
			printer->Density_Data_File("tmp_d-.txt",File_name_MINUScube,"FiniteDifferences_cube");
			//printer->Extrac_N_Write_Integral_DAT("tmp_f-.txt", File_name_MINUScube,"FiniteDifferences_cube");
		}if(fileMinus!="0."){
			Density(Realname, filePlus,File_name_PLUScube,"_dif","+",printer);
			printer->Density_Data_File("tmp_d+.txt",File_name_PLUScube,"FiniteDifferences_cube");
			//printer->Extrac_N_Write_Integral_DAT("tmp_f+.txt", File_name_PLUScube,"FiniteDifferences_cube");
		}if(fileMinus!="0." && filePlus!="0."){
			Density(Realname, fileMinus,File_name_ZEROcube,"_dif","0rad",printer);
			printer->Density_Data_File("tmp_d0rad.txt",File_name_ZEROcube,"FiniteDifferences_cube");
			//printer->Extrac_N_Write_Integral_DAT("tmp_f0.txt", File_name_ZEROcube,"FiniteDifferences_cube");
		}
	}
	//delete tmp files
	printer->DeleteFile("tmp_f+.txt");
	printer->DeleteFile("tmp_f-.txt");
	printer->DeleteFile("tmp_f0.txt");
	printer->DeleteFile("tmp_d+.txt");
	printer->DeleteFile("tmp_d-.txt");
	printer->DeleteFile("tmp_d0rad.txt");
	printer->DeleteFile("tmp.txt");
	//LINUX
	Directory("FiniteDifferences_cube");
	Directory("FiniteDifferences_txt");
	//move files to folders.
	rename(string(File_name_MINUScube+".cube").c_str(),string("FiniteDifferences_cube/f-_dif-"+Realname+".cube").c_str());
	rename(string(File_name_PLUScube+".cube").c_str(),string("FiniteDifferences_cube/f+_dif-"+Realname+".cube").c_str());
	rename(string(File_name_ZEROcube+".cube").c_str(),string("FiniteDifferences_cube/f0rad_dif-"+Realname+".cube").c_str());

	rename(string(File_name_MINUScube+"-VMD.vmd").c_str(),string("FiniteDifferences_cube/f-_dif-"+Realname+"-VMD.vmd").c_str());
	rename(string(File_name_PLUScube+"-VMD.vmd").c_str(),string("FiniteDifferences_cube/f+_dif-"+Realname+"-VMD.vmd").c_str());
	rename(string(File_name_ZEROcube+"-VMD.vmd").c_str(),string("FiniteDifferences_cube/f0rad_dif-"+Realname+"-VMD.vmd").c_str());

	rename(string(File_name_MINUScube+"-integral.txt").c_str(),string("FiniteDifferences_txt/f-_dif-"+Realname+"-integral.txt").c_str());
	rename(string(File_name_PLUScube+"-integral.txt").c_str(),string("FiniteDifferences_txt/f+_dif-"+Realname+"-integral.txt").c_str());
	rename(string(File_name_ZEROcube+"-integral.txt").c_str(),string("FiniteDifferences_txt/f0rad_dif-"+Realname+"-integral.txt").c_str());


	rename(string(File_name_MINUScube+"-rho.cube").c_str(),string("FiniteDifferences_cube/f-_dif-"+Realname+"-rho.cube").c_str());
	rename(string(File_name_PLUScube+"-rho.cube").c_str(),string("FiniteDifferences_cube/f+_dif-"+Realname+"-rho.cube").c_str());
	rename(string(File_name_ZEROcube+"-rho.cube").c_str(),string("FiniteDifferences_cube/f0rad_dif-"+Realname+"-rho.cube").c_str());

	rename(string(File_name_MINUScube+"-rho-VMD.vmd").c_str(),string("FiniteDifferences_cube/f-_dif-"+Realname+"-rho-VMD.vmd").c_str());
	rename(string(File_name_PLUScube+"-rho-VMD.vmd").c_str(),string("FiniteDifferences_cube/f+_dif-"+Realname+"-rho-VMD.vmd").c_str());
	rename(string(File_name_ZEROcube+"-rho-VMD.vmd").c_str(),string("FiniteDifferences_cube/f0rad_dif-"+Realname+"-rho-VMD.vmd").c_str());

	rename(string(File_name_MINUScube+"-rho.txt").c_str(),string("FiniteDifferences_txt/f-_dif-"+Realname+"-rho.txt").c_str());
	rename(string(File_name_PLUScube+"-rho.txt").c_str(),string("FiniteDifferences_txt/f+_dif-"+Realname+"-rho.txt").c_str());
	rename(string(File_name_ZEROcube+"-rho.txt").c_str(),string("FiniteDifferences_txt/f0rad_dif-"+Realname+"-rho.txt").c_str());
}
void approach::Get_Dual(string cubeFPlus,string cubeFMinus,string fileName){
	fileManager *printer=new fileManager();
	string order,dual="Dual-"+fileName+".cube",dualAbs="Dual-"+fileName+"_abs.cube";
	//.cibe pipeline f-
	order="13\n11\n4\n"+cubeFMinus+"\n0\n"+dual+"\nx\n";	
	printer->Multiwfn_Pipeline(order,cubeFPlus.c_str(), "tmpDual.txt");
	order="13\n11\n13\n0\n"+dualAbs+"\nx\n";
	printer->Multiwfn_Pipeline(order,dual.c_str(), "tmpDual.txt");
	order="17\n1\n2\n2\n-1\n"+dual+"\n-3\nx\n";
	printer->Multiwfn_Pipeline(order,dualAbs.c_str(), "tmpDualABS.txt");
	printer->Extrac_N_Write_Integral_DAT("tmpDualABS.txt", "Dual-"+fileName+"_abs","Koopmans_cube");
	printer->DeleteFile("tmpDual.txt");
	printer->DeleteFile("tmpDualABS.txt");
	Directory("Dual_Information");
	rename(string(dual).c_str(),string("Dual_Information/"+dual).c_str());
	rename(string(dualAbs).c_str(),string("Dual_Information/"+dualAbs).c_str());
	rename(string("Dual-"+fileName+"_abs-VMD.vmd").c_str(),string("Dual_Information/Dual-"+fileName+"_abs-VMD.vmd").c_str());
	rename(string("Dual-"+fileName+"_abs-integral.txt").c_str(),string("Dual_Information/Dual-"+fileName+"_abs-integral.txt").c_str());
}

void approach::GetEnergy(ifstream &file, int type){
	string data,sub;
	while(file){
		getline(file,data);
		if(type==0){ //fch file
			size_t f1=data.find("Total Energy");	//secret word to get energy
			if(f1!=string::npos){	
				stringstream iss(data.substr(data.find("R")+1));		//energy is after the R
				iss >> sub;
				break;
			}
			/*if(flag1==1){//we found the position with optimation energy
				stringstream iss(data);
				while(iss){
					iss >> sub;	//we make a array with all the energy data					
					info[contador]=sub;
					contador++;
					if(contador==fch_pos){
						flag1=0;
						break;
					}
				}
				contador--;
			}	
			if(f1!=string::npos){//line found
				stringstream iss(data.substr(data.find("=")+1));	//found the size of the array
				iss >> sub;
				fch_pos= atoi(sub.c_str());
				flag1=1;
				info=new string[fch_pos];
			}*/
		}/*else{//wfn file
			size_t f1=data.find("HF ENERGY");
			if(f1!=string::npos){
				stringstream iss(data.substr(data.find("=")+1)); //info in last line
				iss >> TotalEnergy;
			}
		}*/
	}
	file.close();
	/*
	if(type==0){	//fch file
		for (int i = 0; i < fch_pos; ++i){
			if(i==fch_pos-2){//-1 for array order, -1 for position in file.
				TotalEnergy=strtod(info[i].c_str(), 0);
			}
		}
	}*/
	TotalEnergy=strtod(sub.c_str(),0);	//string to double
	cout.precision(10);
	cout<<"electronic energy:\t"<<TotalEnergy<<endl;
}
void approach::GetHOMO_LUMO(ifstream &file, int type){
	int stop=0;	//this way we read the file until HOMO and LUMO data have been find.
	string sub, data;
	while(file && stop!=2){
		getline(file, data);		//read line by line
		if(type==0){//fch file
			size_t f1=data.find("Number of alpha electrons");		//HOMO's secret password
			size_t f2=data.find("Number of beta electrons");		//LUMO's secret password
			if(f1 != string::npos){	//HOMO
				stringstream iss(data.substr(data.find("I")+1));
				iss >> sub;	//get value
				stop++;
				BossOne[0]=sub;
			}
			if(f2 != string::npos){	//LUMO
				stringstream iss(data.substr(data.find("I")+1));
				iss >> sub;
				stop++;
				BossOne[1]=sub;
			}
		}/*else{	//wfn file
			size_t f1=data.find("MOL ORBITALS");
			if(f1!=string::npos){
				stringstream iss(data.substr(data.find("=")+1));
				iss >> sub;
				iss >> BossOne[0];
				break;
			}
		}*/
	}
	file.close();
	int B1= atoi(BossOne[0].c_str());		//string data to int data.
	int B2= atoi(BossOne[1].c_str());
	ostringstream ss;
	if(B1<B2){	//this get the correct HOMO position
		BossOne[0]=BossOne[1];
		ss << (B2+1);
		BossOne[1]=ss.str();
	}else{
		ss << (B1+1);
		BossOne[1]=ss.str();
	}
}
void approach::GetListCharges(ifstream &file, int type, string txt){
	int flag=0;
	count=0;
	string data;
	charges *helper, *late_one=new charges;	//list obligatory variables and code.
	root->next=NULL;
	root->prev=late_one;
	late_one->next=root;
	helper=root;
	while(file){
		getline(file, data);
		if(type==0){//fch file
			size_t f3=data.find("Alpha Orbital Energies");	//secret block start.
			size_t f4=data.find("Alpha MO coefficients");	//block ends.
			if(f3 !=string::npos || flag==1){	//we found the start or we already found it.
				if (f4 !=string::npos){			//if we find the blocks end then we stop reading the file.
					flag=0;
					break;
				}
				if (flag==1){
					stringstream iss(data);
					while (iss){		//we read each word from te line.
						iss >> helper->charge;		//list is build.
						helper->n_charge=strtod(helper->charge.c_str(), 0);
						helper->next=new charges;
						helper=helper->next;
						late_one=late_one->next;
						helper->prev=late_one;
						count++;
					}
					//this code delete the line jump each 5 chargues.
					helper=late_one;
					late_one=late_one->prev;
					delete helper->next; //free memory.
					helper->next=NULL;
				}
				flag =1;
			}
		}/*else{//WFN file
			//HOMO
			size_t f3=data.find("OCC");
			if(f3!=string::npos){
				aux++;
				stringstream iss(data.substr(data.rfind("=")+1));
				string sub;
				iss >> helper->charge;		//list is build.
				helper->n_charge=strtod(helper->charge.c_str(), 0);
				helper->next=new charges;
				helper=helper->next;
				late_one=late_one->next;
				helper->prev=late_one;
			}
		}*/
	}
	file.close();
	/*if(type!=0 && aux!=atoi(BossOne[0].c_str())){	//WFN file verification
		cout << "WFN file incomplete"<<endl; //file have less chargue than HOMO says.
		exit(1);
	}else if(type!=0){	//WFN file LUMO orbitals chargues
		//LUMO;
		string sub;
		ifstream FILE;
		FILE.open(txt.c_str());
		int verificador=0;
		while(FILE){
			getline(FILE, data);
			stringstream iss(data);
			size_t f2=data.find("END");	//block ends.
			iss >> sub;
			if(aux+1==atoi(sub.c_str()) && verificador==0){
				verificador=1;
				iss >> helper->charge;		//list is build.
					helper->n_charge=strtod(helper->charge.c_str(), 0);
					helper->next=new charges;
					helper=helper->next;
					late_one=late_one->next;
					helper->prev=late_one;
			}else if(verificador!=1){
				cout <<"ERROR: txt file doesn't seem to continue wfn file, HOMO position in wfn file is: ";
				cout <<aux <<"; and LUMO position given is: "<< sub<< "; where "<<aux+1<< " should be"<<endl;	
				exit(1);
			}else{
				if(f2==string::npos){
					stringstream iss(data);
					iss >> helper->charge;		//list is build.
					iss >> helper->charge;		//list is build.
					helper->n_charge=strtod(helper->charge.c_str(), 0);
					helper->next=new charges;
					helper=helper->next;
					late_one=late_one->next;
					helper->prev=late_one;
				}else{
					break;
				}
			}
		}
		FILE.close();
	}
	*/
}
void approach::Show_HOMO_LUMO_Info(){	//simple print.
	int HOMO_position_charge= atoi(BossOne[0].c_str());		//string data to int data.
	int LUMO_position_charge= atoi(BossOne[1].c_str());
	charges *helper = root;	//to read the list with the charges.
	int aux=0;
	if ( helper != 0 ) { //Makes sure there is a place to start
	  	while ( helper->next != 0 ) {
	  		aux++;	//list start with 1.
	    	if(aux==HOMO_position_charge){
	    		BossOne[2]=helper->charge;
	    		homo_position=helper;
	    	}
	    	if(aux==LUMO_position_charge){
	    		BossOne[3]=helper->charge;
	    		lumo_position=helper;
	    	}
	    	helper = helper->next;
	  	}
	}
	cout << "HOMO energy:\t\t\t"<<BossOne[2]<<endl;
	cout << "LUMO energy:\t\t\t" <<BossOne[3]<<endl;
}
void approach::ChemicalPotential(){	//chemical potential using frontier orbitals
	double homo = strtod(BossOne[2].c_str(), 0);
	double lumo = strtod(BossOne[3].c_str(), 0);
	cout.precision(12);
	Che_pot=(homo+lumo)/2;
	cout <<"Chemical potential :\t\t"<<Che_pot<<endl;
	GlobalHardness();	//call global hardness (n)
}

void approach::GlobalHardness(){	//n; using frontier orbitals.
	double homo = strtod(BossOne[2].c_str(), 0);
	double lumo = strtod(BossOne[3].c_str(), 0);
	cout.precision(12);
	float Glo_hard=(lumo-homo)/2;
	cout <<"Global Hardness :\t\t"<<Glo_hard<<endl;
	GlobalSoftness(Glo_hard);
	ElectrophilicityIndex(Glo_hard);
	
}
void approach::ElectrophilicityIndex(float hard){
	float w=(Che_pot*Che_pot)/(2*hard);
	cout <<"Electrophilicity Index (w):\t"<<w<<endl;
}
void approach::Electronegativity(){	//-chemical potential.
	float x=-Che_pot;
	cout <<"Electronegativity :\t\t"<<x<<endl;	
}
void approach::GlobalSoftness(float hard){	//1/global hardness
	float s=1/(2*hard);
	cout <<"Global Softness :\t\t"<<s<<endl;	
	Electronegativity();
}
float approach::IonizationPotencial(){	//-HOMO energy
	double homo = strtod(BossOne[2].c_str(), 0);
	float IP=-homo;
	cout <<"Ionization potential :\t\t"<<IP<<endl;
	return IP;
}
float approach::ElectroAffinity(){		//-LUMO energy
	double lumo = strtod(BossOne[3].c_str(), 0);
	float EA=-lumo;
	cout <<"Electroaffinity :\t\t"<<EA<<endl;
	return EA;
}
void approach::GAP(){	//HOMO energy-LUMO energy
	double homo = strtod(BossOne[2].c_str(), 0);
	double lumo = strtod(BossOne[3].c_str(), 0);
	float GAP=fabs(homo-lumo);
	cout <<"GAP :\t\t\t\t"<<GAP<<endl;
}
float approach::ElectrophilicityDonator(float IP, float EA){ //w- formula
	float nom, denom, square, division;
	nom=(3*IP)+EA;
	denom=(IP-EA)*16;
	square=nom*nom;
	division=square/denom;
	cout <<"w- Electron Donator:\t\t"<<division<<endl;
	return division;
}
float approach::ElectrophilicityAceptor(float IP, float EA){ //w+ formula
	float nom, denom, square, division;
	nom=IP+(EA*3);
	denom=(IP-EA)*16;
	square=nom*nom;
	division=square/denom;
	cout <<"w+ Electron Acceptor:\t\t"<<division<<endl;
	return division;
}
void approach::ElectrophilicityNet(float wMinus, float wPlus){
	float w=wPlus+wMinus;
	cout <<"Electrophilicity Net:\t\t"<<w<<endl;		//w+ + w-
}	


float approach::IonizationPotencialFD(float cation, float neutral){	//Finite differences case, I is cation energy-neutral energy
	float IP=cation-neutral;
	cout <<"Ionization potential :\t\t"<<IP<<endl;
	return IP;
}
float approach::ElectroAffinityFD(float neutral, float anion){ //Finite differences case, I is neutral energy - anion energy
	float EA=neutral-anion;
	cout <<"Electroaffinity :\t\t"<<EA<<endl;
	return EA;
}
void approach::ChemicalPotentialFD(float IP, float EA){	////Finite differences case, u= (I+A)/2
	Che_pot=(-1)*(IP+EA)/2;
	cout <<"Chemical potential :\t\t"<<Che_pot<<endl;
}
void approach::GlobalHardnessFD(float IP, float EA){	//Finite differences case, n= (I-A)
	float n=(IP-EA);
	cout <<"Global Hardness :\t\t"<<n<<endl;
	GlobalSoftness(n);
	ElectrophilicityIndex(n);
}
float approach::getTotalEnergy(){	//get total energy.
	return TotalEnergy;
}
