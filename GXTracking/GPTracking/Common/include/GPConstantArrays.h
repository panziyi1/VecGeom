#ifndef GPConstantArrays_HH
#define GPConstantArrays_HH 1

#include "GPTypeDef.h"

// mean ionisation potential in eV
// IonisationPotential[Z-1]

namespace GPConstantArrays {


CONSTTYPE const size_t numElements = 98;
VARTYPE G4double IonisationPotential[numElements] = 
  {   19.2,  41.8,  40. ,  63.7,  76. ,  81. ,  82. ,  95. , 115. , 137. , 
     149. , 156. , 166. , 173. , 173. , 180. , 174. , 188.0, 190. , 191. ,
     216. , 233. , 245. , 257. , 272. , 286. , 297. , 311. , 322. , 330. ,
     334. , 350. , 347. , 348. , 343. , 352. , 363. , 366. , 379. , 393. ,
     417. , 424. , 428. , 441. , 449. , 470. , 470. , 469. , 488. , 488. ,
     487. , 485. , 491. , 482. , 488. , 491. , 501. , 523. , 535. , 546. ,
     560. , 574. , 580. , 591. , 614. , 628. , 650. , 658. , 674. , 684. ,
     694. , 705. , 718. , 727. , 736. , 746. , 757. , 790. , 790. , 800. ,
     810. , 823. , 823. , 830. , 825. , 794. , 827. , 826. , 841. , 847. ,
     878. , 890. , 902. , 921. , 934. , 939. , 952. , 966. };

};
#endif
