#include <iostream>
#include "menu.h"
#include "aprox.h"

using namespace std;

//menu.h for class definition, variables, etc.
menu::menu(){
//Obligatory constructor
}

void menu::Mainmenu(){
	//vector<approach> ap;
	ifstream FILE;	
	string file=nameFile+"."+extension; //name from ArgumentValidation function chain of command.
	approach *app=new approach();	// class call.
	int n;	//this way we know what to do with HOMO, LUMO positions and charge from files.
	if(extension=="wfn"){	//FUTURE, when we read wfn files this makes the differences in equations.
		n=1;				//for now n=0 always, fch forever
	}else{
		n=0;
	}
	float IP, EA;
	switch(opcion){ //opcion is a int vaue definited in menu.h; it define the approach escoged in argv.
		case 1:	//Koopsman
			FILE.open(file.c_str());
			cout << "Koopman's approximation selected\n" <<endl;
			app->GetHOMO_LUMO(FILE,n);
			FILE.open(file.c_str());
			app->GetListCharges(FILE,n,txtFileName+".txt");
			app->Show_HOMO_LUMO_Info();
			FILE.open(file.c_str());
			app->GAP();
			cout <<"Total ";
			app->GetEnergy(FILE,n);
			cout<<"\n...:::Global reactivity indexes:::...\n\n";
			IP =app->IonizationPotencial();
			EA =app->ElectroAffinity();
			app->ChemicalPotential();
			app->ElectrophilicityNet(app->ElectrophilicityDonator(IP, EA),app->ElectrophilicityAceptor(IP, EA));
			app->Koopmans(file, RealName,work_option,dual);
			break;
		case 2:	//Orbital-weighted Local Reactivity Descriptor
			cout <<"Orbital-Weighted approximation selected\n"<<endl;
			FILE.open(file.c_str());
			app->GetHOMO_LUMO(FILE,n);
			FILE.open(file.c_str());
			app->GetListCharges(FILE,n,txtFileName+".txt");
			app->Show_HOMO_LUMO_Info();
			FILE.open(file.c_str());
			app->GAP();
			cout <<"Total ";
			app->GetEnergy(FILE,n);
			cout<<"\n...:::Global reactivity indexes:::...\n\n";
			IP =app->IonizationPotencial();
			EA =app->ElectroAffinity();
			app->ChemicalPotential();
			app->ElectrophilicityNet(app->ElectrophilicityDonator(IP, EA),app->ElectrophilicityAceptor(IP, EA));
			app->OrbitalWeighted(file,gaussian_value,RealName, work_option);
			break;
		case 3://Finite Differences
			cout <<"Finite differences approximation selected\n"<<endl;
			if(extPlus!="" && extMinus!=""){ //is only possible to calculate Chemical Potential, Softness, Hardness, etc if cation and anion file were given.
				FILE.open((FileNeutral+"."+extension).c_str());
				cout <<"Neutral species ";
				app->GetEnergy(FILE,n);
				float ne=app->getTotalEnergy();
				FILE.open((FileMinus+"."+extMinus).c_str());
				cout <<" Anion species ";
				app->GetEnergy(FILE,n);
				float an=app->getTotalEnergy();
				FILE.open((FilePlus+"."+extPlus).c_str());
				cout <<"Cation species ";
				app->GetEnergy(FILE,n);
				float cat= app->getTotalEnergy();
				//#################################################
				if(cat!=ne || an!=ne){	//if cation, neutral and anion files have the same energy, we dont calculate I, chemical, hardness etc.
					cout<<"\n...:::Global reactivity indexes:::...\n\n";
					IP=app->IonizationPotencialFD(cat,ne);
					EA=app->ElectroAffinityFD(ne,an);
					app->ChemicalPotentialFD(IP,EA);
					app->GlobalHardnessFD(IP,EA);
					app->ElectrophilicityNet(app->ElectrophilicityDonator(IP, EA),app->ElectrophilicityAceptor(IP, EA));
				}
			}
			app->Finite_Differences(FileNeutral+"."+extension, FileMinus+"."+extMinus, FilePlus+"."+extPlus, RealName,work_option);
		default:
			break;
	}
	delete app;
}


void menu::FileValidation(char file[]){	//the file type is supported.
	string fullPath=string(file);
	size_t pos_ext=fullPath.rfind(".");
	extension=fullPath.substr(pos_ext+1);	
	nameFile=fullPath.substr(0,pos_ext);
	if(extension=="wfn" || extension=="fch"){
		FileExist();//we make sure file exist.
		if(extension=="wfn"){
			cout << "Not supported at the moment"<<endl;
			exit(1);
			/*cout <<"Please write txt filename with orbitals energy"<<endl;
			char txt[100];
			cin >> txt;
			txtFileName=nameFile;
			FileValidation(txt);//we validate the txt file given.

		*/}
	/*}else if(extension=="txt"){	//for wfn use
		FileExist();	//we make sure file exist.
		extension="wfn";
		string hlp;	
		hlp=txtFileName;			//change of value in variable to open rigth file in Mainmenu().
		txtFileName=nameFile;
		nameFile=hlp;*/
	}else{
		cout <<"Format file not supported, use -help for use rules"<<endl;
		exit(1);
	}
	size_t pos_slash=nameFile.rfind("/");
	if(pos_slash!=string::npos){						//if the file is a route (/home/user/file.fch)
		RealName= nameFile.substr(pos_slash+1);			//we get the name of the file. 
	}else{												//is later used to name the .cube .dat .vmd files.
		RealName=nameFile;
	}
}
void menu::FileExist(){	//if the file dooesn't open then we close the program.
	string file=nameFile+"."+extension;
	ifstream my_file(file.c_str());
	if (!my_file){
 		cout<<"File given doesn't exist"<<endl;
 		exit(1);
	}
}
void menu::GaussianValidation(char gaussian[]){
	string data=string(gaussian);	//the gaussian value given is a number
	istringstream iss(data);
	iss >> gaussian_value;
	if(iss.eof()==false){
		cout << "Gaussian value isn't a number"<<endl;
		exit(1);
	}
}
void menu::ArgumentsValidation(int argc, char *argv[]){ //Change This function and helps, also the .h file with the new variable work_option
	if(argc==1){	//1 argument will be ./Program.
		cout <<"Insufficients arguments\nDisplay this command-line summary\n\t-help|-h|--help"<<endl;
		exit(1);
	}
	if(!strcmp(argv[1],"-k") || (!strcmp(argv[1], "-koopmans"))){
		if(argc<4 || argc>5){					//make sure input is give.
			DisplayUsageK(argv[0]);		//usage.
		}else{
			if(!strcmp(argv[3],"-d")){
				dual=1;
				FileValidation(argv[4]);	//validate input file.
			}else{
				dual=0;
				FileValidation(argv[3]);	//validate input file.
			}
			opcion=1;					//the file is validated, opcion 1 is Koopman.
		}
	//#//#/#/#/#/#/#/#/#/#//#/#/#/#/#/#/#/#/#/#/#/
	}else if(!strcmp(argv[1],"-w") || (!strcmp(argv[1], "-orbital-weighted"))){
		if(argc!=5){					//4 arguments needed, program -b input gaussian
			DisplayUsageOW(argv[0]);
		}else{
			FileValidation(argv[3]);		//validate file
			GaussianValidation(argv[4]);	//gaussian validation.
		}			
		opcion=2;							//file and gaussian validated, opcion=2 is band.
	//#//#/#/#/#/#/#/#/#/#//#/#/#/#/#/#/#/#/#/#/#/
	}else if(!strcmp(argv[1],"-fd") || (!strcmp(argv[1], "-finite-differences"))){
		if(argc!=6){					//5 arguments needed, program -fd neutraFile cationFile anionFile
			DisplayUsageFD(argv[0]);
		}else{
			
			if(!strcmp(argv[4],"0")){	//cation file not given (0), 
				FileValidation(argv[5]);	//anion file validation.
				FileMinus=nameFile;
				FilePlus="0";
				extPlus="";
				extMinus=extension;
			}else if(!strcmp(argv[5],"0")){	//anion file not given
				FileValidation(argv[4]);	//cation file validation.
				FilePlus=nameFile;
				FileMinus="0";
				extPlus=extension;
				extMinus="";
			}else{							//cation and anion files given
				FileValidation(argv[5]);	//anion file validation.
				FileMinus=nameFile;
				extMinus=extension;
				
				FileValidation(argv[4]);	//cation file validation.
				FilePlus=nameFile;
				extPlus=extension;
			}
			FileValidation(argv[3]);	//neutral file validation.
			FileNeutral=nameFile;
			opcion=3;					//opcion 3 is finite differences

		}
	//#//#/#/#/#/#/#/#/#/#//#/#/#/#/#/#/#/#/#/#/#/
	}else if(!strcmp(argv[1],"-h") || (!strcmp(argv[1], "-help") || (!strcmp(argv[1],"--help")))){
		DisplayHelp(argv[0]);		//show help-
		exit(1);
	}else{
		cout <<"HELP\nDisplay this command-line summary\n\t-help|-h|--help\n"<<endl;
		exit(1);
	}
	if(!strcmp(argv[2],"-f")){
		
		work_option=1;
	}else if(!strcmp(argv[2],"-e")){
		
		work_option=2;
	}else if(!strcmp(argv[2],"-all")){
		
		work_option=3;
	}else{
		cout <<"HELP\nDisplay this command-line summary\n\t-help|-h|--help\n"<<endl;
		exit(1);
	}
}
void menu::DisplayHelp(char name[]){
	cout << "NAME"<<endl;
	cout << "	"<<name<<" - Pipeline for topological analysis of the Fukui function"<<endl;
	cout << ""<<endl;
	cout << "USAGE"<<endl;
	cout << "   -Koopman's approximation\n"<<endl;
	cout << "	 $ "<<name<<" -k -[compute] file.fch\n"<<endl; 
	//cout << "	 $ "<<name<<" -k -[compute] -d file.fch\t\tOptional: Calculate Dual\n"<<endl; 
	cout << "   -Finite difference approximation\n"<<endl;
	cout << "	 $ "<<name<<" -fd -[compute] file-neutral.fch file-cation.fch file-anion.fch\n"<<endl;
	cout << "   Note:"<<endl;
	cout << "   -To obtain the nucleophile Fukui function only, put 0 instead of the anion.fch file."<<endl;
	cout << "   -To obtain the electrophile Fukui function only, put 0 instead of the cation.fch file.\n"<<endl;
	cout << "   -Orbital-weighted Local Reactivity Descriptor approximation\n"<<endl;
	cout << "	 $ "<<name<<" -w -[compute] file.fch [width]\n"<<endl; 
	cout << "	               [width]: Width of the gaussian function (Hartree)\n"<<endl;															
	cout << "\n\nCOMPUTE"<<endl;
	cout << "	-To compute only basins"<<endl;	
	cout << "		[compute]= -f\n"<<endl;	
	cout << "	-To compute only rho"<<endl;	
	cout << "		[compute]= -e\n"<<endl;	
	cout << "	-To compute both"<<endl;	
	cout << "		[compute]= -all\n"<<endl;	

	cout << "BUG REPORT"<<endl;
	cout << "	Report "<<name<<" bugs to osvyanezosses@gmail.com & dinostroza11@alumnos.utalca.cl"<<endl;
	cout << ""<<endl;
	cout << ""<<endl;
}


void menu::DisplayUsageK(char name[]){
	cout << "This script must be run with:" <<endl;
	cout <<"\nUsage:\n"<<name<<" -k -[compute] -[dual] input.fch \n";
	exit(1); 

}
void menu::DisplayUsageOW(char name[]){
	cout << "This script must be run with:" <<endl;
	cout <<"\nUsage:\n"<<name<<" -w -[compute] input.fch [width] \n"; 
	exit(1);

}
void menu::DisplayUsageFD(char name[]){
	cout << "This script must be run with:" <<endl;
	cout <<"\nUsage:\n"<<name<<" -fd -[compute] file-neutral.fch file-cation.fch file-anion.fch \n"; 
	exit(1);

}
