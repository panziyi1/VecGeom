#include <TPartIndex.h>
#include <TString.h>
#include <TDatabasePDG.h>
#include <TMath.h>

ClassImp(TPartIndex)

const char* TPartIndex::fPrName[FNPROC]={"Transport","MultScatt","Ionisation","Decay","inElastic",
			   "Elastic","RestCapture","Brehms","PairProd","Annihilation",
			   "CoulombScatt","Photoel","Compton","Conversion","Capture",
					"Killer","Total"};
const Short_t TPartIndex::fPCode[FNPROC]={1091,2010,2002,6201,4121,4111,4151,2003,2004,2005,2001,
					  2012,2013,2014,4131,7403,999};

static const Int_t kPDG[FNPART]={
521,-521,4232,3322,3312,1000010030,15,-15,4222,4212,3212,3112,3222,2212,-211,211,3334,2112,13,-13,4122,3122,310,130,
-321,321,22,11,-11,1000010020,-4232,-3322,-3312,-1000010030,-4222,-4212,-3112,-3222,-2212,-3334,-2112,-4122,-3122,
-1000010020,-1000020040,-1000020030,1000020040,1000020030,50000060,-431,431,-411,411,20113,215,-215,115,421,-511,-531,
-421,-541,-12212,-12112,-2124,-1214,-22212,-22112,-32212,-32112,-2216,-2116,-12216,-12116,-22124,-21214,-42212,-42112,
-32124,-31214,-42124,-41214,-12218,-12118,-52214,-52114,-2128,-1218,-100002210,-100002110,-100012210,-100012110,531,-5,
-4,-1,-1103,-32214,-32224,-31114,-32114,-2122,-2222,-1112,-1212,-12214,-12224,-11114,-12114,-12122,-12222,-11112,
-11212,-2126,-2226,-1116,-1216,-22122,-22222,-21112,-21212,-22214,-22224,-21114,-22114,-12126,-12226,-11116,-11216,
-2218,-2228,-1118,-2118,-2214,-2224,-1114,-2114,511,-100311,-10311,-10313,-20313,-10315,-315,-100315,-317,-100313,
-30313,-313,-311,541,-13122,-3124,-23122,-33122,-13124,-43122,-53122,-3126,-13126,-23124,-3128,-23126,-5122,443,12212,
-12,-14,-16,12112,-5332,-4332,2124,-3,-3101,-3103,-3224,-3114,-3214,-13222,-13112,-13212,-13224,-13114,-13214,-23222,
-23112,-23212,-3226,-3116,-3216,-13226,-13116,-13216,-23224,-23114,-23214,-3228,-3118,-3218,1214,22212,-3212,-5222,
-5112,-5212,22112,32212,-4112,-3303,-3201,-3203,-6,32112,-2,-2101,-2103,-2203,-3314,-3324,-23314,-23324,-13314,-13324,
-33314,-33324,-13316,-13326,2216,2116,-5132,-5232,12216,-4132,10213,-10213,10113,5,4,50000052,1,1103,32214,32224,31114,
32114,2122,2222,1112,1212,12214,12224,11114,12114,12122,12222,11112,11212,2126,2226,1116,1216,22122,22222,21112,21212,
22214,22224,21114,22114,12126,12226,11116,11216,2218,2228,1118,2118,2214,2224,1114,2114,12116,22124,21214,221,100221,
9020221,100331,10225,10335,331,441,10221,9030221,10331,9000221,9010221,20223,20333,225,9030225,9060225,335,42212,0,21,
10223,10333,100321,-100321,100311,10321,-10321,10311,10323,-10323,10313,20323,-20323,20313,10325,-10325,10315,325,-325,
315,100325,-100325,100315,327,-327,317,100323,-100323,100313,30323,-30323,30313,323,-323,313,42112,32124,311,31214,
42124,41214,13122,3124,23122,33122,13124,43122,53122,3126,13126,23124,3128,23126,5122,12218,12118,52214,52114,12,14,16,
223,100223,30223,2128,227,5332,4332,50000050,333,100333,337,100211,-100211,100111,1218,100002210,111,10215,-10215,
10115,100002110,100213,-100213,100113,30213,-30213,30113,213,-213,113,217,-217,117,3,3101,3103,3224,3114,3214,13222,
13112,13212,13224,13114,13214,23222,23112,23212,3226,3116,3216,13226,13116,13216,23224,23114,23214,3228,3118,3218,
100012210,100012110,553,5222,5112,5212,10211,-10211,4112,3303,3201,3203,6,10111,9000211,-9000211,2,2101,2103,2203,3314,
3324,23314,23324,13314,13324,33314,33324,13316,13326,9000111,20213,5132,5232,-20213,4132};

TPartIndex* TPartIndex::fgPartIndex=0;

//___________________________________________________________________
TPartIndex::TPartIndex():
   fNPart(0),
   fPDG(0),
   fNpReac(0),
   fNEbins(0),
   fEilDelta(0),
   fEGrid(0),
   fDBPdg(TDatabasePDG::Instance())
{ 
}

//___________________________________________________________________
TPartIndex::~TPartIndex() {
   delete [] fPDG;
   delete [] fEGrid;
   delete fDBPdg;
   fgPartIndex=0;
}


//___________________________________________________________________
void TPartIndex::SetEnergyGrid(Double_t emin, Double_t emax, Int_t nbins) {
   fNEbins = nbins;
   fEilDelta = (fNEbins-1)/TMath::Log(emax/emin);
   delete [] fEGrid;
   fEGrid = new Double_t[fNEbins];
   Double_t en=emin;
   Double_t edelta=TMath::Exp(1/fEilDelta);
   for(Int_t i=0; i<fNEbins; ++i) {
      fEGrid[i]=en;
      en*=edelta;
   }
}

//___________________________________________________________________
Int_t TPartIndex::ProcIndex(Int_t proccode) const {
   Short_t ip=fNProc;
   while(ip--) if(fPCode[ip]==proccode) break;
   return ip;
}

//___________________________________________________________________
const Char_t *TPartIndex::ProcName(Int_t proc) const {
   if(proc<0 || proc>=fNProc) return "Unknown";
   return fPrName[proc];
}

//___________________________________________________________________
void TPartIndex::SetPartTable(const Int_t *PDG, Int_t np) {
   fNPart = np;
   fPDG = new Int_t[fNPart];
   for(Int_t i=0; i<fNPart; ++i) 
      fPDG[i]=PDG[i];
}

//______________________________________________________________________________
Int_t TPartIndex::PDG(const Char_t* pname) const {Int_t nr=fNPart;
   while(nr--) if(!strcmp(pname,fDBPdg->GetParticle(fPDG[nr])->GetName())) return fPDG[nr];
   return -12345678;
}

//______________________________________________________________________________
void TPartIndex::Print(Option_t *option) const
{
   Char_t line[120];
   TString opt = option;
   opt.ToLower();
   if(opt.Contains("particles")) {
      printf("Available particles:\n");
      memset(line,0,120);
      for(Int_t i=0; i<fNPart; ++i) {
	 const char *name = fDBPdg->GetParticle(fPDG[i])->GetName();
	 if(strlen(line)+strlen(name)+1>119) {
	    printf("%s\n",line);
	    memset(line,0,120);
	 }
	 strcat(line,name);
	 strcat(line," ");
      }
      if(strlen(line)) printf("%s\n",line);
   }
   if(opt.Contains("reactions")) {
      printf("Available reactions:\n");
      memset(line,0,120);
      strcat(line,"Total ");
      for(Int_t i=0; i<fNProc; ++i) {
	 if(strlen(line)+strlen(fPrName[i])+1>119) {
	    printf("%s\n",line);
	    memset(line,0,120);
	 }
	 strcat(line,fPrName[i]);
	 strcat(line," ");
      }
      if(strlen(line)) printf("%s\n",line);
   }
}

//______________________________________________________________________________
void TPartIndex::Streamer(TBuffer &R__b)
{
   // Stream an object of class TPartIndex.

   if (R__b.IsReading()) {
      delete fDBPdg;
      fDBPdg = 0;
      R__b.ReadClassBuffer(TPartIndex::Class(),this);
      fgPartIndex = this;
   } else {
      R__b.WriteClassBuffer(TPartIndex::Class(),this);
   }
}




