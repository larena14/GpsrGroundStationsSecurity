import inet.common.INETDefs;
import inet.common.TlvOptions;
import inet.common.geometry.Geometry;
import inet.common.packet.chunk.Chunk;
import inet.networklayer.common.L3Address;

cplusplus {{
#include "inet/routing/gpsr/GpsrDefs.h"
}}

namespace inet;


class StationNotice extends FieldsChunk
{
    L3Address address;
    Coord position;
    bool deregister;
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
