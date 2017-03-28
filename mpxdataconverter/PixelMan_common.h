
#pragma once
// common data types
//
typedef char i8;
typedef unsigned char u8;
typedef short i16;
typedef unsigned short u16;
typedef long i32;
typedef unsigned long u32;
typedef int BOOL;
typedef unsigned char byte;
typedef unsigned short DACTYPE;         // type for DAC value
//typedef INT_PTR INTPTR;                 // integral type, also safe for storing pointer (32/64 bit platform)

// variables types identifiers
typedef enum _Data_Types
{
    TYPE_BOOL   = 0,        // C bool value (int)
    TYPE_CHAR   = 1,        // signed char
    TYPE_UCHAR  = 2,        // unsigned char
    TYPE_BYTE   = 3,        // byte (unsigned char)
    TYPE_I16    = 4,        // signed short
    TYPE_U16    = 5,        // unsigned short
    TYPE_I32    = 6,        // int
    TYPE_U32    = 7,        // unsigned int
    TYPE_FLOAT  = 8,        // float
    TYPE_DOUBLE = 9,        // double
    TYPE_STRING = 10,       // zero terminated string
    TYPE_LAST   = 11       // border
} Data_Types;

// structure describing pixel configuration for Mpx
typedef struct _PixelCfg
{
    byte maskBit: 1;        // mask bit (1 bit, low (0) is ACTIVE)
    byte testBit: 1;        // test bit (1 bit, low (0) is ACTIVE)
    byte lowTh: 3;          // low threshold (3 bits, low (0) is ACTIVE)
    byte highTh: 3;         // high threshold (3 bits, low (0) is ACTIVE)
} PixelCfg;

// structure describing pixel configuration for TimePix
typedef struct
{
    byte maskBit: 1;        // mask bit (1 bit, low (0) is ACTIVE)
    byte testBit: 1;        // test bit (1 bit, low (0) is ACTIVE)
    byte thlAdj: 4;         // threshold adjustment
    byte mode: 2;           // pixel mode = {p0, p1}; mode: 0 - medipix, 1 - TOT, 2 - Timepix 1-hit, 3 - Timepix
} TpxPixCfg;


// structure describing "custom" item
typedef struct
{
    Data_Types type;        // variable type
    u32 count;              // array size
    u32 flags;              // combination of flags (e.g. for HwInfoItem MPX_HWINFO_* flags)
    const char *name;       // name of item
    const char *descr;      // item description
    void *data;             // pointer to variable of described type
} ItemInfo;

// structure describing one variable in HW
typedef ItemInfo HwInfoItem;

#define ITEMINFO_CANCHANGE      0x0002      // value can be changed

#define MPX_HWINFO_CFGSAVE      0x0001      // value will be stored in config file
#define MPX_HWINFO_CANCHANGE    0x0002      // value can be changed
#define MPX_HWINFO_MPXFRAME     0x0004      // value be stored for each frame

// macros for type conversions
#define LOWU16(l)           ((u16)((u32)(l) & 0xffff))
#define HIU16(l)            ((u16)((u32)(l) >> 16))
#define MAKEU32(a, b)       ((u32)(((u16)(a)) | (((u32)((u16)(b))) << 16)))

#define MPXCTRL_MAX_MSGLENGTH   4096        // maximum length of error message for some mpxCtrl functions
#define MPX_MAX_MSGLENGTH       4096        // maximum length of error message
#define MPX_MAX_PATH            512         // maximum length of path for file handling
#define MPX_MAX_CHBID           64          // max length of chipboard ID
#define MPX_MAX_IFACENAME       64          // max length of ifaceName name
#define MATRIX_SIZE             65536       // chip matrix size
#define FALSE                   0           // C bool values
#define TRUE                    1

// order of the Medipix DACs in array DACTYPE[] for mpxCtrlSetDACs
typedef enum _DACS_ORDER
{
    DELAYN = 0,
    DISC = 1,
    PREAMP = 2,
    SETDISC = 3,
    THS = 4,
    IKRUM = 5,
    ABUFFER = 6,
    VTHH = 7,
    VTHL = 8,
    VFBK = 9,
    VGND = 10,
    BIASLVDSTX = 11,
    REFLVDSTX = 12,
    IKRUMHALF = 13
} DACS_ORDER;

// order of the Medipix MXR DACs in array DACTYPE[] for mpxCtrlSetDACs
typedef enum _DACS_ORDER_MXR
{
    MXR_IKRUM = 0,
    MXR_DISC = 1,
    MXR_PREAMP = 2,
    MXR_BUFFA = 3,
    MXR_BUFFB = 4,
    MXR_DELAYN = 5,
    MXR_THLFINE = 6,
    MXR_THLCOARSE = 7,
    MXR_THHFINE = 8,
    MXR_THHCOARSE = 9,   
    MXR_FBK = 10,
    MXR_GND = 11,
    MXR_THS = 12,
    MXR_BIASLVDS = 13,
    MXR_REFLVDS = 14
} DACS_ORDER_MXR;

// order of the Medipix TPX DACs in array DACTYPE[] for mpxCtrlSetDACs
typedef enum _DACS_ORDER_TPX
{
    TPX_IKRUM = 0,
    TPX_DISC = 1,
    TPX_PREAMP = 2,
    TPX_BUFFA = 3,
    TPX_BUFFB = 4,
    TPX_HIST = 5,
    TPX_THLFINE = 6,
    TPX_THLCOARSE = 7,
    TPX_VCAS = 8,   
    TPX_FBK = 9,
    TPX_GND = 10,
    TPX_THS = 11,
    TPX_BIASLVDS = 12,
    TPX_REFLVDS = 13
} DACS_ORDER_TPX;

#define MPX_ORIG                    1   // original medipix2.x
#define MPX_MXR                     2   // medipix mxr 2.0
#define MPX_TPX                     3   // timepix

#define TRIGGER_ACQSTART            1   // start trigger (mpxCtrlTrigger(TRIGGER_ACQSTART))
#define TRIGGER_ACQSTOP             2   // stop trigger (mpxCtrlTrigger(TRIGGER_ACQSTOP))

#define ACQMODE_ACQSTART_TIMERSTOP          0x001    // acq. is started immediately (after startAcq call), acq. is stopped by timer (HW or PC) - "common acq"
#define ACQMODE_ACQSTART_HWTRIGSTOP         0x002    // acq. is started immediately (after startAcq call), acq. is stopped by trigger from HW library (e.g. ext shutter monitor)
#define ACQMODE_ACQSTART_SWTRIGSTOP         0x004    // acq. is started immediately (after startAcq call), acq. is stopped by mpxCtrlTrigger(TRIGGER_ACQSTOP) trigger to HW library (e.g. ext shutter monitor)
#define ACQMODE_HWTRIGSTART_TIMERSTOP       0x010    // acq. is started by trigger from HW library, stopped by timer (HW or PC)
#define ACQMODE_HWTRIGSTART_HWTRIGSTOP      0x020    // acq. is started by trigger from HW library, stopped by trigger from HW library
#define ACQMODE_SWTRIGSTART_TIMERSTOP       0x040    // acq. is started by mpxCtrlTrigger(TRIGGER_ACQSTART), stopped by timer  (HW or PC)
#define ACQMODE_SWTRIGSTART_SWTRIGSTOP      0x080    // acq. is started by mpxCtrlTrigger(TRIGGER_ACQSTART), stopped by mpxCtrlTrigger(TRIGGER_ACQSTOP)
#define ACQMODE_EXTSHUTTER                  0x100    // shutter is controlled externally during acq. (e.g. independently on SW)

typedef struct _DevInfo
{
    int pixCount;                       // total number of pixels
    int rowLen;                         // length of row in pixels (e.g 256 for single chip, 512 for quad);
    int numberOfChips;                  // number of chips
    int numberOfRows;                   // number of rows in which chips are aligned (e.g. quad has 4 chips, which are in 2 rows)
    int mpxType;                        // medipix type - MPX_ORIG, MPX_MXR, MPX_TPX
    char chipboardID[MPX_MAX_CHBID];    // id of chip/chipboard
    const char *ifaceName;              // name of interface

    u32 suppAcqModes;                   // supported acq. mode bitwise combinations of ACQMODE_xxx flags
    BOOL suppCallback;                  // callback is supported (acq. is finished, triggers)

    double clockFreq;                   // clock frequency [MHz] (going to MPX) - for Timepix information
    // hw timer capabilities
    double maxTimer;                    // maximum of hw timer [s]
    double timeElement;                 // timer elemenet (minimum nonzero time) [s]

    // test pulse capabilities
    u32 maxPulseCount;                  // maximum number of pulses that can be sent
    double maxPulseHeight;              // max pulse height [V]
    double maxPulsePeriod;              // max period of pulses [s], length of pulse should be period/2

    // ext DAC capabilities
    double extDacMinV;                  // minimum external DAC voltage
    double extDacMaxV;                  // maximum external DAC voltage
    double extDacStep;                  // ext. DAC step size
} DevInfo;


// format flags for saving to file (MgrSaveFrameType or flags for perform acq. functions)
// flags FSAVE_BINARY/FSAVE_ASCII and flags FSAVE_I16/FSAVE_U32/FSAVE_DOUBLE are mutually exclusive
// if FSAVE_SPARSEX and FSAVE_SPARSEXY is not specified full matrix is saved
#define FSAVE_BINARY        0x0001      // save in binary format
#define FSAVE_ASCII         0x0002      // save in ASCII format
#define FSAVE_APPEND        0x0004      // append data file to existing file if exists
#define FSAVE_I16           0x0010      // save as 16bit integer
#define FSAVE_U32           0x0020      // save as unsigned 32bit integer
#define FSAVE_DOUBLE        0x0040      // save as double
#define FSAVE_NODESCFILE    0x0100      // do not save description file
#define FSAVE_SPARSEXY      0x1000      // save only nonzero position in [x y count] format
#define FSAVE_SPARSEX       0x2000      // save only nonzero position in [x count] format
#define FSAVE_NOFILE        0x8000      // frame will not be saved :)


#define sstrncpy(strDest, strSource, count)\
    strncpy(strDest, strSource, count), (strDest)[count-1] = '\0', (strDest)

static const u32 sizeofType[TYPE_LAST] =
{
    sizeof(BOOL),
    sizeof(char),
    sizeof(unsigned char),
    sizeof(byte),
    sizeof(i16),
    sizeof(u16),
    sizeof(i32),
    sizeof(u32),
    sizeof(float),
    sizeof(double),
    sizeof(char*)
};

static const char * const nameofType[] = 
{
    "bool",
    "char",
    "uchar",
    "byte",
    "i16",
    "u16",
    "i32",
    "u32",
    "float",
    "double",
    "string",
    "uknown"
};
