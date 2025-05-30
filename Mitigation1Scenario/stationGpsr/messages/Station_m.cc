//
// Generated file, do not edit! Created by opp_msgtool 6.0 from stationGpsr/messages/Station.msg.
//

// Disable warnings about unused variables, empty switch stmts, etc:
#ifdef _MSC_VER
#  pragma warning(disable:4101)
#  pragma warning(disable:4065)
#endif

#if defined(__clang__)
#  pragma clang diagnostic ignored "-Wshadow"
#  pragma clang diagnostic ignored "-Wconversion"
#  pragma clang diagnostic ignored "-Wunused-parameter"
#  pragma clang diagnostic ignored "-Wc++98-compat"
#  pragma clang diagnostic ignored "-Wunreachable-code-break"
#  pragma clang diagnostic ignored "-Wold-style-cast"
#elif defined(__GNUC__)
#  pragma GCC diagnostic ignored "-Wshadow"
#  pragma GCC diagnostic ignored "-Wconversion"
#  pragma GCC diagnostic ignored "-Wunused-parameter"
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#  pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#  pragma GCC diagnostic ignored "-Wfloat-conversion"
#endif

#include <iostream>
#include <sstream>
#include <memory>
#include <type_traits>
#include "Station_m.h"

namespace omnetpp {

// Template pack/unpack rules. They are declared *after* a1l type-specific pack functions for multiple reasons.
// They are in the omnetpp namespace, to allow them to be found by argument-dependent lookup via the cCommBuffer argument

// Packing/unpacking an std::vector
template<typename T, typename A>
void doParsimPacking(omnetpp::cCommBuffer *buffer, const std::vector<T,A>& v)
{
    int n = v.size();
    doParsimPacking(buffer, n);
    for (int i = 0; i < n; i++)
        doParsimPacking(buffer, v[i]);
}

template<typename T, typename A>
void doParsimUnpacking(omnetpp::cCommBuffer *buffer, std::vector<T,A>& v)
{
    int n;
    doParsimUnpacking(buffer, n);
    v.resize(n);
    for (int i = 0; i < n; i++)
        doParsimUnpacking(buffer, v[i]);
}

// Packing/unpacking an std::list
template<typename T, typename A>
void doParsimPacking(omnetpp::cCommBuffer *buffer, const std::list<T,A>& l)
{
    doParsimPacking(buffer, (int)l.size());
    for (typename std::list<T,A>::const_iterator it = l.begin(); it != l.end(); ++it)
        doParsimPacking(buffer, (T&)*it);
}

template<typename T, typename A>
void doParsimUnpacking(omnetpp::cCommBuffer *buffer, std::list<T,A>& l)
{
    int n;
    doParsimUnpacking(buffer, n);
    for (int i = 0; i < n; i++) {
        l.push_back(T());
        doParsimUnpacking(buffer, l.back());
    }
}

// Packing/unpacking an std::set
template<typename T, typename Tr, typename A>
void doParsimPacking(omnetpp::cCommBuffer *buffer, const std::set<T,Tr,A>& s)
{
    doParsimPacking(buffer, (int)s.size());
    for (typename std::set<T,Tr,A>::const_iterator it = s.begin(); it != s.end(); ++it)
        doParsimPacking(buffer, *it);
}

template<typename T, typename Tr, typename A>
void doParsimUnpacking(omnetpp::cCommBuffer *buffer, std::set<T,Tr,A>& s)
{
    int n;
    doParsimUnpacking(buffer, n);
    for (int i = 0; i < n; i++) {
        T x;
        doParsimUnpacking(buffer, x);
        s.insert(x);
    }
}

// Packing/unpacking an std::map
template<typename K, typename V, typename Tr, typename A>
void doParsimPacking(omnetpp::cCommBuffer *buffer, const std::map<K,V,Tr,A>& m)
{
    doParsimPacking(buffer, (int)m.size());
    for (typename std::map<K,V,Tr,A>::const_iterator it = m.begin(); it != m.end(); ++it) {
        doParsimPacking(buffer, it->first);
        doParsimPacking(buffer, it->second);
    }
}

template<typename K, typename V, typename Tr, typename A>
void doParsimUnpacking(omnetpp::cCommBuffer *buffer, std::map<K,V,Tr,A>& m)
{
    int n;
    doParsimUnpacking(buffer, n);
    for (int i = 0; i < n; i++) {
        K k; V v;
        doParsimUnpacking(buffer, k);
        doParsimUnpacking(buffer, v);
        m[k] = v;
    }
}

// Default pack/unpack function for arrays
template<typename T>
void doParsimArrayPacking(omnetpp::cCommBuffer *b, const T *t, int n)
{
    for (int i = 0; i < n; i++)
        doParsimPacking(b, t[i]);
}

template<typename T>
void doParsimArrayUnpacking(omnetpp::cCommBuffer *b, T *t, int n)
{
    for (int i = 0; i < n; i++)
        doParsimUnpacking(b, t[i]);
}

// Default rule to prevent compiler from choosing base class' doParsimPacking() function
template<typename T>
void doParsimPacking(omnetpp::cCommBuffer *, const T& t)
{
    throw omnetpp::cRuntimeError("Parsim error: No doParsimPacking() function for type %s", omnetpp::opp_typename(typeid(t)));
}

template<typename T>
void doParsimUnpacking(omnetpp::cCommBuffer *, T& t)
{
    throw omnetpp::cRuntimeError("Parsim error: No doParsimUnpacking() function for type %s", omnetpp::opp_typename(typeid(t)));
}

}  // namespace omnetpp

namespace inet {

Register_Class(StationNotice)

StationNotice::StationNotice() : ::inet::FieldsChunk()
{
}

StationNotice::StationNotice(const StationNotice& other) : ::inet::FieldsChunk(other)
{
    copy(other);
}

StationNotice::~StationNotice()
{
}

StationNotice& StationNotice::operator=(const StationNotice& other)
{
    if (this == &other) return *this;
    ::inet::FieldsChunk::operator=(other);
    copy(other);
    return *this;
}

void StationNotice::copy(const StationNotice& other)
{
    this->source = other.source;
    this->sourceModuleName = other.sourceModuleName;
    this->position = other.position;
    this->deregister = other.deregister;
}

void StationNotice::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::FieldsChunk::parsimPack(b);
    doParsimPacking(b,this->source);
    doParsimPacking(b,this->sourceModuleName);
    doParsimPacking(b,this->position);
    doParsimPacking(b,this->deregister);
}

void StationNotice::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::FieldsChunk::parsimUnpack(b);
    doParsimUnpacking(b,this->source);
    doParsimUnpacking(b,this->sourceModuleName);
    doParsimUnpacking(b,this->position);
    doParsimUnpacking(b,this->deregister);
}

const L3Address& StationNotice::getSource() const
{
    return this->source;
}

void StationNotice::setSource(const L3Address& source)
{
    handleChange();
    this->source = source;
}

const char * StationNotice::getSourceModuleName() const
{
    return this->sourceModuleName.c_str();
}

void StationNotice::setSourceModuleName(const char * sourceModuleName)
{
    handleChange();
    this->sourceModuleName = sourceModuleName;
}

const Coord& StationNotice::getPosition() const
{
    return this->position;
}

void StationNotice::setPosition(const Coord& position)
{
    handleChange();
    this->position = position;
}

bool StationNotice::getDeregister() const
{
    return this->deregister;
}

void StationNotice::setDeregister(bool deregister)
{
    handleChange();
    this->deregister = deregister;
}

class StationNoticeDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_source,
        FIELD_sourceModuleName,
        FIELD_position,
        FIELD_deregister,
    };
  public:
    StationNoticeDescriptor();
    virtual ~StationNoticeDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyName) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyName) const override;
    virtual int getFieldArraySize(omnetpp::any_ptr object, int field) const override;
    virtual void setFieldArraySize(omnetpp::any_ptr object, int field, int size) const override;

    virtual const char *getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const override;
    virtual std::string getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const override;
    virtual omnetpp::cValue getFieldValue(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual omnetpp::any_ptr getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const override;
};

Register_ClassDescriptor(StationNoticeDescriptor)

StationNoticeDescriptor::StationNoticeDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::StationNotice)), "inet::FieldsChunk")
{
    propertyNames = nullptr;
}

StationNoticeDescriptor::~StationNoticeDescriptor()
{
    delete[] propertyNames;
}

bool StationNoticeDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<StationNotice *>(obj)!=nullptr;
}

const char **StationNoticeDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *StationNoticeDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int StationNoticeDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 4+base->getFieldCount() : 4;
}

unsigned int StationNoticeDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        0,    // FIELD_source
        FD_ISEDITABLE,    // FIELD_sourceModuleName
        FD_ISCOMPOUND,    // FIELD_position
        FD_ISEDITABLE,    // FIELD_deregister
    };
    return (field >= 0 && field < 4) ? fieldTypeFlags[field] : 0;
}

const char *StationNoticeDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "source",
        "sourceModuleName",
        "position",
        "deregister",
    };
    return (field >= 0 && field < 4) ? fieldNames[field] : nullptr;
}

int StationNoticeDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "source") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "sourceModuleName") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "position") == 0) return baseIndex + 2;
    if (strcmp(fieldName, "deregister") == 0) return baseIndex + 3;
    return base ? base->findField(fieldName) : -1;
}

const char *StationNoticeDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "inet::L3Address",    // FIELD_source
        "string",    // FIELD_sourceModuleName
        "inet::Coord",    // FIELD_position
        "bool",    // FIELD_deregister
    };
    return (field >= 0 && field < 4) ? fieldTypeStrings[field] : nullptr;
}

const char **StationNoticeDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldPropertyNames(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

const char *StationNoticeDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

int StationNoticeDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    StationNotice *pp = omnetpp::fromAnyPtr<StationNotice>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void StationNoticeDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    StationNotice *pp = omnetpp::fromAnyPtr<StationNotice>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'StationNotice'", field);
    }
}

const char *StationNoticeDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    StationNotice *pp = omnetpp::fromAnyPtr<StationNotice>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string StationNoticeDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    StationNotice *pp = omnetpp::fromAnyPtr<StationNotice>(object); (void)pp;
    switch (field) {
        case FIELD_source: return pp->getSource().str();
        case FIELD_sourceModuleName: return oppstring2string(pp->getSourceModuleName());
        case FIELD_position: return pp->getPosition().str();
        case FIELD_deregister: return bool2string(pp->getDeregister());
        default: return "";
    }
}

void StationNoticeDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    StationNotice *pp = omnetpp::fromAnyPtr<StationNotice>(object); (void)pp;
    switch (field) {
        case FIELD_sourceModuleName: pp->setSourceModuleName((value)); break;
        case FIELD_deregister: pp->setDeregister(string2bool(value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'StationNotice'", field);
    }
}

omnetpp::cValue StationNoticeDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    StationNotice *pp = omnetpp::fromAnyPtr<StationNotice>(object); (void)pp;
    switch (field) {
        case FIELD_source: return omnetpp::toAnyPtr(&pp->getSource()); break;
        case FIELD_sourceModuleName: return pp->getSourceModuleName();
        case FIELD_position: return omnetpp::toAnyPtr(&pp->getPosition()); break;
        case FIELD_deregister: return pp->getDeregister();
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'StationNotice' as cValue -- field index out of range?", field);
    }
}

void StationNoticeDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    StationNotice *pp = omnetpp::fromAnyPtr<StationNotice>(object); (void)pp;
    switch (field) {
        case FIELD_sourceModuleName: pp->setSourceModuleName(value.stringValue()); break;
        case FIELD_deregister: pp->setDeregister(value.boolValue()); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'StationNotice'", field);
    }
}

const char *StationNoticeDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructName(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        case FIELD_position: return omnetpp::opp_typename(typeid(Coord));
        default: return nullptr;
    };
}

omnetpp::any_ptr StationNoticeDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    StationNotice *pp = omnetpp::fromAnyPtr<StationNotice>(object); (void)pp;
    switch (field) {
        case FIELD_source: return omnetpp::toAnyPtr(&pp->getSource()); break;
        case FIELD_position: return omnetpp::toAnyPtr(&pp->getPosition()); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void StationNoticeDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    StationNotice *pp = omnetpp::fromAnyPtr<StationNotice>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'StationNotice'", field);
    }
}

Register_Class(StationNoticeResponse)

StationNoticeResponse::StationNoticeResponse() : ::inet::FieldsChunk()
{
}

StationNoticeResponse::StationNoticeResponse(const StationNoticeResponse& other) : ::inet::FieldsChunk(other)
{
    copy(other);
}

StationNoticeResponse::~StationNoticeResponse()
{
}

StationNoticeResponse& StationNoticeResponse::operator=(const StationNoticeResponse& other)
{
    if (this == &other) return *this;
    ::inet::FieldsChunk::operator=(other);
    copy(other);
    return *this;
}

void StationNoticeResponse::copy(const StationNoticeResponse& other)
{
    this->c = other.c;
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    this->nonce = other.nonce;
}

void StationNoticeResponse::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::FieldsChunk::parsimPack(b);
    doParsimPacking(b,this->c);
    doParsimPacking(b,this->x);
    doParsimPacking(b,this->y);
    doParsimPacking(b,this->z);
    doParsimPacking(b,this->nonce);
}

void StationNoticeResponse::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::FieldsChunk::parsimUnpack(b);
    doParsimUnpacking(b,this->c);
    doParsimUnpacking(b,this->x);
    doParsimUnpacking(b,this->y);
    doParsimUnpacking(b,this->z);
    doParsimUnpacking(b,this->nonce);
}

const char * StationNoticeResponse::getC() const
{
    return this->c.c_str();
}

void StationNoticeResponse::setC(const char * c)
{
    handleChange();
    this->c = c;
}

double StationNoticeResponse::getX() const
{
    return this->x;
}

void StationNoticeResponse::setX(double x)
{
    handleChange();
    this->x = x;
}

double StationNoticeResponse::getY() const
{
    return this->y;
}

void StationNoticeResponse::setY(double y)
{
    handleChange();
    this->y = y;
}

double StationNoticeResponse::getZ() const
{
    return this->z;
}

void StationNoticeResponse::setZ(double z)
{
    handleChange();
    this->z = z;
}

uint64_t StationNoticeResponse::getNonce() const
{
    return this->nonce;
}

void StationNoticeResponse::setNonce(uint64_t nonce)
{
    handleChange();
    this->nonce = nonce;
}

class StationNoticeResponseDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_c,
        FIELD_x,
        FIELD_y,
        FIELD_z,
        FIELD_nonce,
    };
  public:
    StationNoticeResponseDescriptor();
    virtual ~StationNoticeResponseDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyName) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyName) const override;
    virtual int getFieldArraySize(omnetpp::any_ptr object, int field) const override;
    virtual void setFieldArraySize(omnetpp::any_ptr object, int field, int size) const override;

    virtual const char *getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const override;
    virtual std::string getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const override;
    virtual omnetpp::cValue getFieldValue(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual omnetpp::any_ptr getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const override;
};

Register_ClassDescriptor(StationNoticeResponseDescriptor)

StationNoticeResponseDescriptor::StationNoticeResponseDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::StationNoticeResponse)), "inet::FieldsChunk")
{
    propertyNames = nullptr;
}

StationNoticeResponseDescriptor::~StationNoticeResponseDescriptor()
{
    delete[] propertyNames;
}

bool StationNoticeResponseDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<StationNoticeResponse *>(obj)!=nullptr;
}

const char **StationNoticeResponseDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *StationNoticeResponseDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int StationNoticeResponseDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 5+base->getFieldCount() : 5;
}

unsigned int StationNoticeResponseDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISEDITABLE,    // FIELD_c
        FD_ISEDITABLE,    // FIELD_x
        FD_ISEDITABLE,    // FIELD_y
        FD_ISEDITABLE,    // FIELD_z
        FD_ISEDITABLE,    // FIELD_nonce
    };
    return (field >= 0 && field < 5) ? fieldTypeFlags[field] : 0;
}

const char *StationNoticeResponseDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "c",
        "x",
        "y",
        "z",
        "nonce",
    };
    return (field >= 0 && field < 5) ? fieldNames[field] : nullptr;
}

int StationNoticeResponseDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "c") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "x") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "y") == 0) return baseIndex + 2;
    if (strcmp(fieldName, "z") == 0) return baseIndex + 3;
    if (strcmp(fieldName, "nonce") == 0) return baseIndex + 4;
    return base ? base->findField(fieldName) : -1;
}

const char *StationNoticeResponseDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "string",    // FIELD_c
        "double",    // FIELD_x
        "double",    // FIELD_y
        "double",    // FIELD_z
        "uint64_t",    // FIELD_nonce
    };
    return (field >= 0 && field < 5) ? fieldTypeStrings[field] : nullptr;
}

const char **StationNoticeResponseDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldPropertyNames(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

const char *StationNoticeResponseDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

int StationNoticeResponseDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    StationNoticeResponse *pp = omnetpp::fromAnyPtr<StationNoticeResponse>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void StationNoticeResponseDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    StationNoticeResponse *pp = omnetpp::fromAnyPtr<StationNoticeResponse>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'StationNoticeResponse'", field);
    }
}

const char *StationNoticeResponseDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    StationNoticeResponse *pp = omnetpp::fromAnyPtr<StationNoticeResponse>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string StationNoticeResponseDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    StationNoticeResponse *pp = omnetpp::fromAnyPtr<StationNoticeResponse>(object); (void)pp;
    switch (field) {
        case FIELD_c: return oppstring2string(pp->getC());
        case FIELD_x: return double2string(pp->getX());
        case FIELD_y: return double2string(pp->getY());
        case FIELD_z: return double2string(pp->getZ());
        case FIELD_nonce: return uint642string(pp->getNonce());
        default: return "";
    }
}

void StationNoticeResponseDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    StationNoticeResponse *pp = omnetpp::fromAnyPtr<StationNoticeResponse>(object); (void)pp;
    switch (field) {
        case FIELD_c: pp->setC((value)); break;
        case FIELD_x: pp->setX(string2double(value)); break;
        case FIELD_y: pp->setY(string2double(value)); break;
        case FIELD_z: pp->setZ(string2double(value)); break;
        case FIELD_nonce: pp->setNonce(string2uint64(value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'StationNoticeResponse'", field);
    }
}

omnetpp::cValue StationNoticeResponseDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    StationNoticeResponse *pp = omnetpp::fromAnyPtr<StationNoticeResponse>(object); (void)pp;
    switch (field) {
        case FIELD_c: return pp->getC();
        case FIELD_x: return pp->getX();
        case FIELD_y: return pp->getY();
        case FIELD_z: return pp->getZ();
        case FIELD_nonce: return (omnetpp::intval_t)(pp->getNonce());
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'StationNoticeResponse' as cValue -- field index out of range?", field);
    }
}

void StationNoticeResponseDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    StationNoticeResponse *pp = omnetpp::fromAnyPtr<StationNoticeResponse>(object); (void)pp;
    switch (field) {
        case FIELD_c: pp->setC(value.stringValue()); break;
        case FIELD_x: pp->setX(value.doubleValue()); break;
        case FIELD_y: pp->setY(value.doubleValue()); break;
        case FIELD_z: pp->setZ(value.doubleValue()); break;
        case FIELD_nonce: pp->setNonce(omnetpp::checked_int_cast<uint64_t>(value.intValue())); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'StationNoticeResponse'", field);
    }
}

const char *StationNoticeResponseDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructName(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    };
}

omnetpp::any_ptr StationNoticeResponseDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    StationNoticeResponse *pp = omnetpp::fromAnyPtr<StationNoticeResponse>(object); (void)pp;
    switch (field) {
        default: return omnetpp::any_ptr(nullptr);
    }
}

void StationNoticeResponseDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    StationNoticeResponse *pp = omnetpp::fromAnyPtr<StationNoticeResponse>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'StationNoticeResponse'", field);
    }
}

Register_Class(PositionRequest)

PositionRequest::PositionRequest() : ::inet::FieldsChunk()
{
}

PositionRequest::PositionRequest(const PositionRequest& other) : ::inet::FieldsChunk(other)
{
    copy(other);
}

PositionRequest::~PositionRequest()
{
}

PositionRequest& PositionRequest::operator=(const PositionRequest& other)
{
    if (this == &other) return *this;
    ::inet::FieldsChunk::operator=(other);
    copy(other);
    return *this;
}

void PositionRequest::copy(const PositionRequest& other)
{
    this->source = other.source;
    this->sourceModuleName = other.sourceModuleName;
    this->address = other.address;
}

void PositionRequest::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::FieldsChunk::parsimPack(b);
    doParsimPacking(b,this->source);
    doParsimPacking(b,this->sourceModuleName);
    doParsimPacking(b,this->address);
}

void PositionRequest::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::FieldsChunk::parsimUnpack(b);
    doParsimUnpacking(b,this->source);
    doParsimUnpacking(b,this->sourceModuleName);
    doParsimUnpacking(b,this->address);
}

const L3Address& PositionRequest::getSource() const
{
    return this->source;
}

void PositionRequest::setSource(const L3Address& source)
{
    handleChange();
    this->source = source;
}

const char * PositionRequest::getSourceModuleName() const
{
    return this->sourceModuleName.c_str();
}

void PositionRequest::setSourceModuleName(const char * sourceModuleName)
{
    handleChange();
    this->sourceModuleName = sourceModuleName;
}

const L3Address& PositionRequest::getAddress() const
{
    return this->address;
}

void PositionRequest::setAddress(const L3Address& address)
{
    handleChange();
    this->address = address;
}

class PositionRequestDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_source,
        FIELD_sourceModuleName,
        FIELD_address,
    };
  public:
    PositionRequestDescriptor();
    virtual ~PositionRequestDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyName) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyName) const override;
    virtual int getFieldArraySize(omnetpp::any_ptr object, int field) const override;
    virtual void setFieldArraySize(omnetpp::any_ptr object, int field, int size) const override;

    virtual const char *getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const override;
    virtual std::string getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const override;
    virtual omnetpp::cValue getFieldValue(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual omnetpp::any_ptr getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const override;
};

Register_ClassDescriptor(PositionRequestDescriptor)

PositionRequestDescriptor::PositionRequestDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::PositionRequest)), "inet::FieldsChunk")
{
    propertyNames = nullptr;
}

PositionRequestDescriptor::~PositionRequestDescriptor()
{
    delete[] propertyNames;
}

bool PositionRequestDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<PositionRequest *>(obj)!=nullptr;
}

const char **PositionRequestDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *PositionRequestDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int PositionRequestDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 3+base->getFieldCount() : 3;
}

unsigned int PositionRequestDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        0,    // FIELD_source
        FD_ISEDITABLE,    // FIELD_sourceModuleName
        0,    // FIELD_address
    };
    return (field >= 0 && field < 3) ? fieldTypeFlags[field] : 0;
}

const char *PositionRequestDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "source",
        "sourceModuleName",
        "address",
    };
    return (field >= 0 && field < 3) ? fieldNames[field] : nullptr;
}

int PositionRequestDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "source") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "sourceModuleName") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "address") == 0) return baseIndex + 2;
    return base ? base->findField(fieldName) : -1;
}

const char *PositionRequestDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "inet::L3Address",    // FIELD_source
        "string",    // FIELD_sourceModuleName
        "inet::L3Address",    // FIELD_address
    };
    return (field >= 0 && field < 3) ? fieldTypeStrings[field] : nullptr;
}

const char **PositionRequestDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldPropertyNames(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

const char *PositionRequestDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

int PositionRequestDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    PositionRequest *pp = omnetpp::fromAnyPtr<PositionRequest>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void PositionRequestDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    PositionRequest *pp = omnetpp::fromAnyPtr<PositionRequest>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'PositionRequest'", field);
    }
}

const char *PositionRequestDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    PositionRequest *pp = omnetpp::fromAnyPtr<PositionRequest>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string PositionRequestDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    PositionRequest *pp = omnetpp::fromAnyPtr<PositionRequest>(object); (void)pp;
    switch (field) {
        case FIELD_source: return pp->getSource().str();
        case FIELD_sourceModuleName: return oppstring2string(pp->getSourceModuleName());
        case FIELD_address: return pp->getAddress().str();
        default: return "";
    }
}

void PositionRequestDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    PositionRequest *pp = omnetpp::fromAnyPtr<PositionRequest>(object); (void)pp;
    switch (field) {
        case FIELD_sourceModuleName: pp->setSourceModuleName((value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'PositionRequest'", field);
    }
}

omnetpp::cValue PositionRequestDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    PositionRequest *pp = omnetpp::fromAnyPtr<PositionRequest>(object); (void)pp;
    switch (field) {
        case FIELD_source: return omnetpp::toAnyPtr(&pp->getSource()); break;
        case FIELD_sourceModuleName: return pp->getSourceModuleName();
        case FIELD_address: return omnetpp::toAnyPtr(&pp->getAddress()); break;
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'PositionRequest' as cValue -- field index out of range?", field);
    }
}

void PositionRequestDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    PositionRequest *pp = omnetpp::fromAnyPtr<PositionRequest>(object); (void)pp;
    switch (field) {
        case FIELD_sourceModuleName: pp->setSourceModuleName(value.stringValue()); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'PositionRequest'", field);
    }
}

const char *PositionRequestDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructName(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    };
}

omnetpp::any_ptr PositionRequestDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    PositionRequest *pp = omnetpp::fromAnyPtr<PositionRequest>(object); (void)pp;
    switch (field) {
        case FIELD_source: return omnetpp::toAnyPtr(&pp->getSource()); break;
        case FIELD_address: return omnetpp::toAnyPtr(&pp->getAddress()); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void PositionRequestDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    PositionRequest *pp = omnetpp::fromAnyPtr<PositionRequest>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'PositionRequest'", field);
    }
}

Register_Class(PositionResponse)

PositionResponse::PositionResponse() : ::inet::FieldsChunk()
{
}

PositionResponse::PositionResponse(const PositionResponse& other) : ::inet::FieldsChunk(other)
{
    copy(other);
}

PositionResponse::~PositionResponse()
{
}

PositionResponse& PositionResponse::operator=(const PositionResponse& other)
{
    if (this == &other) return *this;
    ::inet::FieldsChunk::operator=(other);
    copy(other);
    return *this;
}

void PositionResponse::copy(const PositionResponse& other)
{
    this->setted = other.setted;
    this->address = other.address;
    this->position = other.position;
    this->time = other.time;
}

void PositionResponse::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::FieldsChunk::parsimPack(b);
    doParsimPacking(b,this->setted);
    doParsimPacking(b,this->address);
    doParsimPacking(b,this->position);
    doParsimPacking(b,this->time);
}

void PositionResponse::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::FieldsChunk::parsimUnpack(b);
    doParsimUnpacking(b,this->setted);
    doParsimUnpacking(b,this->address);
    doParsimUnpacking(b,this->position);
    doParsimUnpacking(b,this->time);
}

bool PositionResponse::getSetted() const
{
    return this->setted;
}

void PositionResponse::setSetted(bool setted)
{
    handleChange();
    this->setted = setted;
}

const L3Address& PositionResponse::getAddress() const
{
    return this->address;
}

void PositionResponse::setAddress(const L3Address& address)
{
    handleChange();
    this->address = address;
}

const Coord& PositionResponse::getPosition() const
{
    return this->position;
}

void PositionResponse::setPosition(const Coord& position)
{
    handleChange();
    this->position = position;
}

::omnetpp::simtime_t PositionResponse::getTime() const
{
    return this->time;
}

void PositionResponse::setTime(::omnetpp::simtime_t time)
{
    handleChange();
    this->time = time;
}

class PositionResponseDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_setted,
        FIELD_address,
        FIELD_position,
        FIELD_time,
    };
  public:
    PositionResponseDescriptor();
    virtual ~PositionResponseDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyName) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyName) const override;
    virtual int getFieldArraySize(omnetpp::any_ptr object, int field) const override;
    virtual void setFieldArraySize(omnetpp::any_ptr object, int field, int size) const override;

    virtual const char *getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const override;
    virtual std::string getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const override;
    virtual omnetpp::cValue getFieldValue(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual omnetpp::any_ptr getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const override;
};

Register_ClassDescriptor(PositionResponseDescriptor)

PositionResponseDescriptor::PositionResponseDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::PositionResponse)), "inet::FieldsChunk")
{
    propertyNames = nullptr;
}

PositionResponseDescriptor::~PositionResponseDescriptor()
{
    delete[] propertyNames;
}

bool PositionResponseDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<PositionResponse *>(obj)!=nullptr;
}

const char **PositionResponseDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *PositionResponseDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int PositionResponseDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 4+base->getFieldCount() : 4;
}

unsigned int PositionResponseDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISEDITABLE,    // FIELD_setted
        0,    // FIELD_address
        FD_ISCOMPOUND,    // FIELD_position
        FD_ISEDITABLE,    // FIELD_time
    };
    return (field >= 0 && field < 4) ? fieldTypeFlags[field] : 0;
}

const char *PositionResponseDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "setted",
        "address",
        "position",
        "time",
    };
    return (field >= 0 && field < 4) ? fieldNames[field] : nullptr;
}

int PositionResponseDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "setted") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "address") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "position") == 0) return baseIndex + 2;
    if (strcmp(fieldName, "time") == 0) return baseIndex + 3;
    return base ? base->findField(fieldName) : -1;
}

const char *PositionResponseDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "bool",    // FIELD_setted
        "inet::L3Address",    // FIELD_address
        "inet::Coord",    // FIELD_position
        "omnetpp::simtime_t",    // FIELD_time
    };
    return (field >= 0 && field < 4) ? fieldTypeStrings[field] : nullptr;
}

const char **PositionResponseDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldPropertyNames(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

const char *PositionResponseDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

int PositionResponseDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    PositionResponse *pp = omnetpp::fromAnyPtr<PositionResponse>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void PositionResponseDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    PositionResponse *pp = omnetpp::fromAnyPtr<PositionResponse>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'PositionResponse'", field);
    }
}

const char *PositionResponseDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    PositionResponse *pp = omnetpp::fromAnyPtr<PositionResponse>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string PositionResponseDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    PositionResponse *pp = omnetpp::fromAnyPtr<PositionResponse>(object); (void)pp;
    switch (field) {
        case FIELD_setted: return bool2string(pp->getSetted());
        case FIELD_address: return pp->getAddress().str();
        case FIELD_position: return pp->getPosition().str();
        case FIELD_time: return simtime2string(pp->getTime());
        default: return "";
    }
}

void PositionResponseDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    PositionResponse *pp = omnetpp::fromAnyPtr<PositionResponse>(object); (void)pp;
    switch (field) {
        case FIELD_setted: pp->setSetted(string2bool(value)); break;
        case FIELD_time: pp->setTime(string2simtime(value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'PositionResponse'", field);
    }
}

omnetpp::cValue PositionResponseDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    PositionResponse *pp = omnetpp::fromAnyPtr<PositionResponse>(object); (void)pp;
    switch (field) {
        case FIELD_setted: return pp->getSetted();
        case FIELD_address: return omnetpp::toAnyPtr(&pp->getAddress()); break;
        case FIELD_position: return omnetpp::toAnyPtr(&pp->getPosition()); break;
        case FIELD_time: return pp->getTime().dbl();
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'PositionResponse' as cValue -- field index out of range?", field);
    }
}

void PositionResponseDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    PositionResponse *pp = omnetpp::fromAnyPtr<PositionResponse>(object); (void)pp;
    switch (field) {
        case FIELD_setted: pp->setSetted(value.boolValue()); break;
        case FIELD_time: pp->setTime(value.doubleValue()); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'PositionResponse'", field);
    }
}

const char *PositionResponseDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructName(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        case FIELD_position: return omnetpp::opp_typename(typeid(Coord));
        default: return nullptr;
    };
}

omnetpp::any_ptr PositionResponseDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    PositionResponse *pp = omnetpp::fromAnyPtr<PositionResponse>(object); (void)pp;
    switch (field) {
        case FIELD_address: return omnetpp::toAnyPtr(&pp->getAddress()); break;
        case FIELD_position: return omnetpp::toAnyPtr(&pp->getPosition()); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void PositionResponseDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    PositionResponse *pp = omnetpp::fromAnyPtr<PositionResponse>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'PositionResponse'", field);
    }
}

Register_Class(S2SPositionRequest)

S2SPositionRequest::S2SPositionRequest() : ::inet::FieldsChunk()
{
}

S2SPositionRequest::S2SPositionRequest(const S2SPositionRequest& other) : ::inet::FieldsChunk(other)
{
    copy(other);
}

S2SPositionRequest::~S2SPositionRequest()
{
}

S2SPositionRequest& S2SPositionRequest::operator=(const S2SPositionRequest& other)
{
    if (this == &other) return *this;
    ::inet::FieldsChunk::operator=(other);
    copy(other);
    return *this;
}

void S2SPositionRequest::copy(const S2SPositionRequest& other)
{
    this->applicant = other.applicant;
    this->applicantModuleName = other.applicantModuleName;
    this->source = other.source;
    this->sourceModuleName = other.sourceModuleName;
    this->address = other.address;
}

void S2SPositionRequest::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::FieldsChunk::parsimPack(b);
    doParsimPacking(b,this->applicant);
    doParsimPacking(b,this->applicantModuleName);
    doParsimPacking(b,this->source);
    doParsimPacking(b,this->sourceModuleName);
    doParsimPacking(b,this->address);
}

void S2SPositionRequest::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::FieldsChunk::parsimUnpack(b);
    doParsimUnpacking(b,this->applicant);
    doParsimUnpacking(b,this->applicantModuleName);
    doParsimUnpacking(b,this->source);
    doParsimUnpacking(b,this->sourceModuleName);
    doParsimUnpacking(b,this->address);
}

const L3Address& S2SPositionRequest::getApplicant() const
{
    return this->applicant;
}

void S2SPositionRequest::setApplicant(const L3Address& applicant)
{
    handleChange();
    this->applicant = applicant;
}

const char * S2SPositionRequest::getApplicantModuleName() const
{
    return this->applicantModuleName.c_str();
}

void S2SPositionRequest::setApplicantModuleName(const char * applicantModuleName)
{
    handleChange();
    this->applicantModuleName = applicantModuleName;
}

const L3Address& S2SPositionRequest::getSource() const
{
    return this->source;
}

void S2SPositionRequest::setSource(const L3Address& source)
{
    handleChange();
    this->source = source;
}

const char * S2SPositionRequest::getSourceModuleName() const
{
    return this->sourceModuleName.c_str();
}

void S2SPositionRequest::setSourceModuleName(const char * sourceModuleName)
{
    handleChange();
    this->sourceModuleName = sourceModuleName;
}

const L3Address& S2SPositionRequest::getAddress() const
{
    return this->address;
}

void S2SPositionRequest::setAddress(const L3Address& address)
{
    handleChange();
    this->address = address;
}

class S2SPositionRequestDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_applicant,
        FIELD_applicantModuleName,
        FIELD_source,
        FIELD_sourceModuleName,
        FIELD_address,
    };
  public:
    S2SPositionRequestDescriptor();
    virtual ~S2SPositionRequestDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyName) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyName) const override;
    virtual int getFieldArraySize(omnetpp::any_ptr object, int field) const override;
    virtual void setFieldArraySize(omnetpp::any_ptr object, int field, int size) const override;

    virtual const char *getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const override;
    virtual std::string getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const override;
    virtual omnetpp::cValue getFieldValue(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual omnetpp::any_ptr getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const override;
};

Register_ClassDescriptor(S2SPositionRequestDescriptor)

S2SPositionRequestDescriptor::S2SPositionRequestDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::S2SPositionRequest)), "inet::FieldsChunk")
{
    propertyNames = nullptr;
}

S2SPositionRequestDescriptor::~S2SPositionRequestDescriptor()
{
    delete[] propertyNames;
}

bool S2SPositionRequestDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<S2SPositionRequest *>(obj)!=nullptr;
}

const char **S2SPositionRequestDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *S2SPositionRequestDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int S2SPositionRequestDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 5+base->getFieldCount() : 5;
}

unsigned int S2SPositionRequestDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        0,    // FIELD_applicant
        FD_ISEDITABLE,    // FIELD_applicantModuleName
        0,    // FIELD_source
        FD_ISEDITABLE,    // FIELD_sourceModuleName
        0,    // FIELD_address
    };
    return (field >= 0 && field < 5) ? fieldTypeFlags[field] : 0;
}

const char *S2SPositionRequestDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "applicant",
        "applicantModuleName",
        "source",
        "sourceModuleName",
        "address",
    };
    return (field >= 0 && field < 5) ? fieldNames[field] : nullptr;
}

int S2SPositionRequestDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "applicant") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "applicantModuleName") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "source") == 0) return baseIndex + 2;
    if (strcmp(fieldName, "sourceModuleName") == 0) return baseIndex + 3;
    if (strcmp(fieldName, "address") == 0) return baseIndex + 4;
    return base ? base->findField(fieldName) : -1;
}

const char *S2SPositionRequestDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "inet::L3Address",    // FIELD_applicant
        "string",    // FIELD_applicantModuleName
        "inet::L3Address",    // FIELD_source
        "string",    // FIELD_sourceModuleName
        "inet::L3Address",    // FIELD_address
    };
    return (field >= 0 && field < 5) ? fieldTypeStrings[field] : nullptr;
}

const char **S2SPositionRequestDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldPropertyNames(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

const char *S2SPositionRequestDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

int S2SPositionRequestDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    S2SPositionRequest *pp = omnetpp::fromAnyPtr<S2SPositionRequest>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void S2SPositionRequestDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    S2SPositionRequest *pp = omnetpp::fromAnyPtr<S2SPositionRequest>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'S2SPositionRequest'", field);
    }
}

const char *S2SPositionRequestDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    S2SPositionRequest *pp = omnetpp::fromAnyPtr<S2SPositionRequest>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string S2SPositionRequestDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    S2SPositionRequest *pp = omnetpp::fromAnyPtr<S2SPositionRequest>(object); (void)pp;
    switch (field) {
        case FIELD_applicant: return pp->getApplicant().str();
        case FIELD_applicantModuleName: return oppstring2string(pp->getApplicantModuleName());
        case FIELD_source: return pp->getSource().str();
        case FIELD_sourceModuleName: return oppstring2string(pp->getSourceModuleName());
        case FIELD_address: return pp->getAddress().str();
        default: return "";
    }
}

void S2SPositionRequestDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    S2SPositionRequest *pp = omnetpp::fromAnyPtr<S2SPositionRequest>(object); (void)pp;
    switch (field) {
        case FIELD_applicantModuleName: pp->setApplicantModuleName((value)); break;
        case FIELD_sourceModuleName: pp->setSourceModuleName((value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'S2SPositionRequest'", field);
    }
}

omnetpp::cValue S2SPositionRequestDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    S2SPositionRequest *pp = omnetpp::fromAnyPtr<S2SPositionRequest>(object); (void)pp;
    switch (field) {
        case FIELD_applicant: return omnetpp::toAnyPtr(&pp->getApplicant()); break;
        case FIELD_applicantModuleName: return pp->getApplicantModuleName();
        case FIELD_source: return omnetpp::toAnyPtr(&pp->getSource()); break;
        case FIELD_sourceModuleName: return pp->getSourceModuleName();
        case FIELD_address: return omnetpp::toAnyPtr(&pp->getAddress()); break;
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'S2SPositionRequest' as cValue -- field index out of range?", field);
    }
}

void S2SPositionRequestDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    S2SPositionRequest *pp = omnetpp::fromAnyPtr<S2SPositionRequest>(object); (void)pp;
    switch (field) {
        case FIELD_applicantModuleName: pp->setApplicantModuleName(value.stringValue()); break;
        case FIELD_sourceModuleName: pp->setSourceModuleName(value.stringValue()); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'S2SPositionRequest'", field);
    }
}

const char *S2SPositionRequestDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructName(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    };
}

omnetpp::any_ptr S2SPositionRequestDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    S2SPositionRequest *pp = omnetpp::fromAnyPtr<S2SPositionRequest>(object); (void)pp;
    switch (field) {
        case FIELD_applicant: return omnetpp::toAnyPtr(&pp->getApplicant()); break;
        case FIELD_source: return omnetpp::toAnyPtr(&pp->getSource()); break;
        case FIELD_address: return omnetpp::toAnyPtr(&pp->getAddress()); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void S2SPositionRequestDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    S2SPositionRequest *pp = omnetpp::fromAnyPtr<S2SPositionRequest>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'S2SPositionRequest'", field);
    }
}

Register_Class(S2SPositionResponse)

S2SPositionResponse::S2SPositionResponse() : ::inet::FieldsChunk()
{
}

S2SPositionResponse::S2SPositionResponse(const S2SPositionResponse& other) : ::inet::FieldsChunk(other)
{
    copy(other);
}

S2SPositionResponse::~S2SPositionResponse()
{
}

S2SPositionResponse& S2SPositionResponse::operator=(const S2SPositionResponse& other)
{
    if (this == &other) return *this;
    ::inet::FieldsChunk::operator=(other);
    copy(other);
    return *this;
}

void S2SPositionResponse::copy(const S2SPositionResponse& other)
{
    this->applicant = other.applicant;
    this->applicantModuleName = other.applicantModuleName;
    this->setted = other.setted;
    this->address = other.address;
    this->position = other.position;
    this->time = other.time;
}

void S2SPositionResponse::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::FieldsChunk::parsimPack(b);
    doParsimPacking(b,this->applicant);
    doParsimPacking(b,this->applicantModuleName);
    doParsimPacking(b,this->setted);
    doParsimPacking(b,this->address);
    doParsimPacking(b,this->position);
    doParsimPacking(b,this->time);
}

void S2SPositionResponse::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::FieldsChunk::parsimUnpack(b);
    doParsimUnpacking(b,this->applicant);
    doParsimUnpacking(b,this->applicantModuleName);
    doParsimUnpacking(b,this->setted);
    doParsimUnpacking(b,this->address);
    doParsimUnpacking(b,this->position);
    doParsimUnpacking(b,this->time);
}

const L3Address& S2SPositionResponse::getApplicant() const
{
    return this->applicant;
}

void S2SPositionResponse::setApplicant(const L3Address& applicant)
{
    handleChange();
    this->applicant = applicant;
}

const char * S2SPositionResponse::getApplicantModuleName() const
{
    return this->applicantModuleName.c_str();
}

void S2SPositionResponse::setApplicantModuleName(const char * applicantModuleName)
{
    handleChange();
    this->applicantModuleName = applicantModuleName;
}

bool S2SPositionResponse::getSetted() const
{
    return this->setted;
}

void S2SPositionResponse::setSetted(bool setted)
{
    handleChange();
    this->setted = setted;
}

const L3Address& S2SPositionResponse::getAddress() const
{
    return this->address;
}

void S2SPositionResponse::setAddress(const L3Address& address)
{
    handleChange();
    this->address = address;
}

const Coord& S2SPositionResponse::getPosition() const
{
    return this->position;
}

void S2SPositionResponse::setPosition(const Coord& position)
{
    handleChange();
    this->position = position;
}

::omnetpp::simtime_t S2SPositionResponse::getTime() const
{
    return this->time;
}

void S2SPositionResponse::setTime(::omnetpp::simtime_t time)
{
    handleChange();
    this->time = time;
}

class S2SPositionResponseDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_applicant,
        FIELD_applicantModuleName,
        FIELD_setted,
        FIELD_address,
        FIELD_position,
        FIELD_time,
    };
  public:
    S2SPositionResponseDescriptor();
    virtual ~S2SPositionResponseDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyName) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyName) const override;
    virtual int getFieldArraySize(omnetpp::any_ptr object, int field) const override;
    virtual void setFieldArraySize(omnetpp::any_ptr object, int field, int size) const override;

    virtual const char *getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const override;
    virtual std::string getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const override;
    virtual omnetpp::cValue getFieldValue(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual omnetpp::any_ptr getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const override;
};

Register_ClassDescriptor(S2SPositionResponseDescriptor)

S2SPositionResponseDescriptor::S2SPositionResponseDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::S2SPositionResponse)), "inet::FieldsChunk")
{
    propertyNames = nullptr;
}

S2SPositionResponseDescriptor::~S2SPositionResponseDescriptor()
{
    delete[] propertyNames;
}

bool S2SPositionResponseDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<S2SPositionResponse *>(obj)!=nullptr;
}

const char **S2SPositionResponseDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *S2SPositionResponseDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int S2SPositionResponseDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 6+base->getFieldCount() : 6;
}

unsigned int S2SPositionResponseDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        0,    // FIELD_applicant
        FD_ISEDITABLE,    // FIELD_applicantModuleName
        FD_ISEDITABLE,    // FIELD_setted
        0,    // FIELD_address
        FD_ISCOMPOUND,    // FIELD_position
        FD_ISEDITABLE,    // FIELD_time
    };
    return (field >= 0 && field < 6) ? fieldTypeFlags[field] : 0;
}

const char *S2SPositionResponseDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "applicant",
        "applicantModuleName",
        "setted",
        "address",
        "position",
        "time",
    };
    return (field >= 0 && field < 6) ? fieldNames[field] : nullptr;
}

int S2SPositionResponseDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "applicant") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "applicantModuleName") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "setted") == 0) return baseIndex + 2;
    if (strcmp(fieldName, "address") == 0) return baseIndex + 3;
    if (strcmp(fieldName, "position") == 0) return baseIndex + 4;
    if (strcmp(fieldName, "time") == 0) return baseIndex + 5;
    return base ? base->findField(fieldName) : -1;
}

const char *S2SPositionResponseDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "inet::L3Address",    // FIELD_applicant
        "string",    // FIELD_applicantModuleName
        "bool",    // FIELD_setted
        "inet::L3Address",    // FIELD_address
        "inet::Coord",    // FIELD_position
        "omnetpp::simtime_t",    // FIELD_time
    };
    return (field >= 0 && field < 6) ? fieldTypeStrings[field] : nullptr;
}

const char **S2SPositionResponseDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldPropertyNames(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

const char *S2SPositionResponseDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

int S2SPositionResponseDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    S2SPositionResponse *pp = omnetpp::fromAnyPtr<S2SPositionResponse>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void S2SPositionResponseDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    S2SPositionResponse *pp = omnetpp::fromAnyPtr<S2SPositionResponse>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'S2SPositionResponse'", field);
    }
}

const char *S2SPositionResponseDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    S2SPositionResponse *pp = omnetpp::fromAnyPtr<S2SPositionResponse>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string S2SPositionResponseDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    S2SPositionResponse *pp = omnetpp::fromAnyPtr<S2SPositionResponse>(object); (void)pp;
    switch (field) {
        case FIELD_applicant: return pp->getApplicant().str();
        case FIELD_applicantModuleName: return oppstring2string(pp->getApplicantModuleName());
        case FIELD_setted: return bool2string(pp->getSetted());
        case FIELD_address: return pp->getAddress().str();
        case FIELD_position: return pp->getPosition().str();
        case FIELD_time: return simtime2string(pp->getTime());
        default: return "";
    }
}

void S2SPositionResponseDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    S2SPositionResponse *pp = omnetpp::fromAnyPtr<S2SPositionResponse>(object); (void)pp;
    switch (field) {
        case FIELD_applicantModuleName: pp->setApplicantModuleName((value)); break;
        case FIELD_setted: pp->setSetted(string2bool(value)); break;
        case FIELD_time: pp->setTime(string2simtime(value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'S2SPositionResponse'", field);
    }
}

omnetpp::cValue S2SPositionResponseDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    S2SPositionResponse *pp = omnetpp::fromAnyPtr<S2SPositionResponse>(object); (void)pp;
    switch (field) {
        case FIELD_applicant: return omnetpp::toAnyPtr(&pp->getApplicant()); break;
        case FIELD_applicantModuleName: return pp->getApplicantModuleName();
        case FIELD_setted: return pp->getSetted();
        case FIELD_address: return omnetpp::toAnyPtr(&pp->getAddress()); break;
        case FIELD_position: return omnetpp::toAnyPtr(&pp->getPosition()); break;
        case FIELD_time: return pp->getTime().dbl();
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'S2SPositionResponse' as cValue -- field index out of range?", field);
    }
}

void S2SPositionResponseDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    S2SPositionResponse *pp = omnetpp::fromAnyPtr<S2SPositionResponse>(object); (void)pp;
    switch (field) {
        case FIELD_applicantModuleName: pp->setApplicantModuleName(value.stringValue()); break;
        case FIELD_setted: pp->setSetted(value.boolValue()); break;
        case FIELD_time: pp->setTime(value.doubleValue()); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'S2SPositionResponse'", field);
    }
}

const char *S2SPositionResponseDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructName(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        case FIELD_position: return omnetpp::opp_typename(typeid(Coord));
        default: return nullptr;
    };
}

omnetpp::any_ptr S2SPositionResponseDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    S2SPositionResponse *pp = omnetpp::fromAnyPtr<S2SPositionResponse>(object); (void)pp;
    switch (field) {
        case FIELD_applicant: return omnetpp::toAnyPtr(&pp->getApplicant()); break;
        case FIELD_address: return omnetpp::toAnyPtr(&pp->getAddress()); break;
        case FIELD_position: return omnetpp::toAnyPtr(&pp->getPosition()); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void S2SPositionResponseDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    S2SPositionResponse *pp = omnetpp::fromAnyPtr<S2SPositionResponse>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'S2SPositionResponse'", field);
    }
}

Register_Class(SimpleMessage)

SimpleMessage::SimpleMessage() : ::inet::FieldsChunk()
{
}

SimpleMessage::SimpleMessage(const SimpleMessage& other) : ::inet::FieldsChunk(other)
{
    copy(other);
}

SimpleMessage::~SimpleMessage()
{
}

SimpleMessage& SimpleMessage::operator=(const SimpleMessage& other)
{
    if (this == &other) return *this;
    ::inet::FieldsChunk::operator=(other);
    copy(other);
    return *this;
}

void SimpleMessage::copy(const SimpleMessage& other)
{
    this->payload = other.payload;
}

void SimpleMessage::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::FieldsChunk::parsimPack(b);
    doParsimPacking(b,this->payload);
}

void SimpleMessage::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::FieldsChunk::parsimUnpack(b);
    doParsimUnpacking(b,this->payload);
}

const char * SimpleMessage::getPayload() const
{
    return this->payload.c_str();
}

void SimpleMessage::setPayload(const char * payload)
{
    handleChange();
    this->payload = payload;
}

class SimpleMessageDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_payload,
    };
  public:
    SimpleMessageDescriptor();
    virtual ~SimpleMessageDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyName) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyName) const override;
    virtual int getFieldArraySize(omnetpp::any_ptr object, int field) const override;
    virtual void setFieldArraySize(omnetpp::any_ptr object, int field, int size) const override;

    virtual const char *getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const override;
    virtual std::string getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const override;
    virtual omnetpp::cValue getFieldValue(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual omnetpp::any_ptr getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const override;
};

Register_ClassDescriptor(SimpleMessageDescriptor)

SimpleMessageDescriptor::SimpleMessageDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::SimpleMessage)), "inet::FieldsChunk")
{
    propertyNames = nullptr;
}

SimpleMessageDescriptor::~SimpleMessageDescriptor()
{
    delete[] propertyNames;
}

bool SimpleMessageDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<SimpleMessage *>(obj)!=nullptr;
}

const char **SimpleMessageDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *SimpleMessageDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int SimpleMessageDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 1+base->getFieldCount() : 1;
}

unsigned int SimpleMessageDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISEDITABLE,    // FIELD_payload
    };
    return (field >= 0 && field < 1) ? fieldTypeFlags[field] : 0;
}

const char *SimpleMessageDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "payload",
    };
    return (field >= 0 && field < 1) ? fieldNames[field] : nullptr;
}

int SimpleMessageDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "payload") == 0) return baseIndex + 0;
    return base ? base->findField(fieldName) : -1;
}

const char *SimpleMessageDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "string",    // FIELD_payload
    };
    return (field >= 0 && field < 1) ? fieldTypeStrings[field] : nullptr;
}

const char **SimpleMessageDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldPropertyNames(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

const char *SimpleMessageDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

int SimpleMessageDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    SimpleMessage *pp = omnetpp::fromAnyPtr<SimpleMessage>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void SimpleMessageDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    SimpleMessage *pp = omnetpp::fromAnyPtr<SimpleMessage>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'SimpleMessage'", field);
    }
}

const char *SimpleMessageDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    SimpleMessage *pp = omnetpp::fromAnyPtr<SimpleMessage>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string SimpleMessageDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    SimpleMessage *pp = omnetpp::fromAnyPtr<SimpleMessage>(object); (void)pp;
    switch (field) {
        case FIELD_payload: return oppstring2string(pp->getPayload());
        default: return "";
    }
}

void SimpleMessageDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    SimpleMessage *pp = omnetpp::fromAnyPtr<SimpleMessage>(object); (void)pp;
    switch (field) {
        case FIELD_payload: pp->setPayload((value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'SimpleMessage'", field);
    }
}

omnetpp::cValue SimpleMessageDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    SimpleMessage *pp = omnetpp::fromAnyPtr<SimpleMessage>(object); (void)pp;
    switch (field) {
        case FIELD_payload: return pp->getPayload();
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'SimpleMessage' as cValue -- field index out of range?", field);
    }
}

void SimpleMessageDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    SimpleMessage *pp = omnetpp::fromAnyPtr<SimpleMessage>(object); (void)pp;
    switch (field) {
        case FIELD_payload: pp->setPayload(value.stringValue()); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'SimpleMessage'", field);
    }
}

const char *SimpleMessageDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructName(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    };
}

omnetpp::any_ptr SimpleMessageDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    SimpleMessage *pp = omnetpp::fromAnyPtr<SimpleMessage>(object); (void)pp;
    switch (field) {
        default: return omnetpp::any_ptr(nullptr);
    }
}

void SimpleMessageDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    SimpleMessage *pp = omnetpp::fromAnyPtr<SimpleMessage>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'SimpleMessage'", field);
    }
}

}  // namespace inet

namespace omnetpp {

}  // namespace omnetpp

