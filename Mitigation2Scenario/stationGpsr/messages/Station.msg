import inet.common.INETDefs;
import inet.common.TlvOptions;
import inet.common.geometry.Geometry;
import inet.common.packet.chunk.Chunk;
import inet.networklayer.common.L3Address;

cplusplus {{

    
#include "inet/routing/gpsr/GpsrDefs.h"
}}

namespace inet;

//MODIFICA VA AGGIUNTO IL CAMPO NONCE E LA LISTA DEGLI INDIRIZZI DEI SUOI VICINI:
class StationNotice extends FieldsChunk
{
    L3Address source;//modifica
    string sourceModuleName; //modifica
    Coord position;
    bool deregister;
    //uint64_t nonce; //modifica
    L3Address addressList[];
    double distanceList[];
    Coord positionList[];
}

class verifyPositionRequest extends FieldsChunk
{
    L3Address source;
    L3Address addressToVerify;
}    
class verifyPositionResponse extends FieldsChunk
{
    L3Address source;
    L3Address AddressToVerify;
    int status; //0 non è presente 1 è presente
}



// MODIFICA AGGIUNTA MESSAGGI RISPOSTA STAZIONE STATIONNOTICE

class StationNoticeResponse extends FieldsChunk
{
    string c; // CIFRATURA
    double x;
    double y;
    double z;
    uint64_t nonce; //NONCE
}

class NeighborVerificationRequest extends FieldsChunk {
    L3Address target; // Indirizzo del nodo da verificare
}

class PositionRequest extends FieldsChunk
{
    L3Address source;
    string sourceModuleName;
    L3Address address;
}

// vedere se eliminare timestamp
class PositionResponse extends FieldsChunk
{
    bool setted;
    L3Address address;
    Coord position;
    simtime_t time; 
}

class S2SPositionRequest extends FieldsChunk
{
    L3Address applicant;
    string applicantModuleName;
    L3Address source;
    string sourceModuleName;
    L3Address address;
}

// vedere se eliminare timestamp
class S2SPositionResponse extends FieldsChunk
{
    L3Address applicant;
    string applicantModuleName;
    bool setted;
    L3Address address;
    Coord position;
    simtime_t time;
}

class SimpleMessage extends FieldsChunk
{
    string payload;    
}
