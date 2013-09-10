#ifndef KRAKENCONFIG_H
#define KRAKENCONFIG_H

#include "src/Kraken/KrakenMap.h"

class KrakenConfig
{
 public:
  KrakenConfig(Kraken * p);

  bool Configure(const string & fileName);

 private:
  KrakenConfig() {
    m_pKraken = NULL;
  }

  Kraken * m_pKraken;

};



#endif //

