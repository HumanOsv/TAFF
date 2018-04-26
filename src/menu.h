#ifndef MENU_H
#define MENU_H


//definicion de la clase validaciones...
class menu {
	std::string nameFile;	//file name given, with all the route.
	std::string txtFileName;	//used in wfn files
	std::string RealName;		//name use to write f+ f- f0rad
	std::string extension;		//fch or wfn
	std::string extMinus;		//finite differences case, anion extension
	std::string extPlus;		//finite differences case, cation extension
	std::string FileNeutral;	//finite differences case, neutral route
	std::string FileMinus;		//finite differences case, anion route
	std::string FilePlus;		//finite differences case, cation route.
	int opcion, work_option,dual;					//opcion 0=koopsman, 1=band 2=fd
	float gaussian_value;		//gausian value.
	//declaracion de todo lo publico...
	public:
	
		menu();//constructor de la clase...
		void Mainmenu();
		void ArgumentsValidation(int, char *argv[]);
		void FileValidation(char[]);
		void FileExist();
		void GaussianValidation(char[]);
		void DisplayUsageK(char[]);
		void DisplayUsageOW(char[]);
		void DisplayUsageFD(char[]);
		void DisplayHelp(char[]);
		
		
};//fin definicion de la clase...
	
#endif

