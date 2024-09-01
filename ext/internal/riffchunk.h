#ifndef INTERNAL_RIFFCHUNK_H
#define INTERNAL_RIFFCHUNK_H

// For detail, see:
// https://web.archive.org/web/20080113195252/http://www.borg.com/~jglatt/tech/wave.htm
typedef uint32_t  ChunkID;
#define MakeFOURCC( ch0, ch1, ch2, ch3 ) \
 ( (uint32_t)(uint8_t)(ch0) | ( (uint32_t)(uint8_t)(ch1) << 8 ) | \
 ( (uint32_t)(uint8_t)(ch2) << 16 ) | ( (uint32_t)(uint8_t)(ch3) << 24 ) )

#define FOURCC_RIFF            MakeFOURCC('R', 'I', 'F', 'F')
#define FOURCC_WAVE            MakeFOURCC('W', 'A', 'V', 'E')
#define ChunkID_Format         MakeFOURCC('f', 'm', 't', ' ')
#define ChunkID_Data           MakeFOURCC('d', 'a', 't', 'a')
#define ChunkID_Cue            MakeFOURCC('c', 'u', 'e', ' ')
#define ChunkID_Playlist       MakeFOURCC('p', 'l', 's', 't')
#define ChunkID_List           MakeFOURCC('l', 'i', 's', 't')
#define ChunkID_AssocDataList  MakeFOURCC('a', 'd', 't', 'l')
#define ChunkID_LABEL          MakeFOURCC('l', 'a', 'b', 'l')
#define ChunkID_NOTE           MakeFOURCC('n', 'o', 't', 'e')
#define ChunkID_LabeledText    MakeFOURCC('l', 't', 'x', 't')
#define ChunkID_Sample         MakeFOURCC('s', 'm', 'p', 'l')
#define ChunkID_Instrument     MakeFOURCC('i', 'n', 's', 't')


/* Format Chunk */
typedef struct {
	ChunkID   chunk_ID; /* 'fmt ' */
	int32_t   chunk_size;

	int16_t   format_tag;
	uint16_t  channels;
	uint32_t  samples_per_sec;
	uint32_t  bytes_per_sec;
	uint16_t  block_size;
	uint16_t  bits_per_sample;

} FormatChunk;

/* Data Chunk */
typedef struct {
	ChunkID   chunk_ID; /* 'data' */
	int32_t   chunk_Size;

	uint8_t   waveform_data[];
} DataChunk;


/* Cue Chunk */
typedef struct {
	int32_t    identifier;
	int32_t    position;
	ChunkID    fcc_chunk;
	int32_t    chunk_start;
	int32_t    block_Start;
	int32_t    sample_offset;
} CuePoint;

typedef struct {
	ChunkID    chunk_ID; /* 'cue ' */
	int32_t    chunk_size;

	int32_t    cue_points;
	CuePoint   points[];
} CueChunk;


/* Playlist Chunk */
typedef struct {
	int32_t    identifier;
	int32_t    length;
	int32_t    repeats;
} Segment;

typedef struct {
	ChunkID    chunk_ID; /* 'plst' */
	int32_t    chunk_size;

	int32_t    dwSegments;
	Segment    segments[];
} PlaylistChunk;


/* Associated Data List */
typedef struct {
	ChunkID    list_ID;  /* 'list' */
	int32_t    chunk_size;  /* includes the Type ID below */
	ChunkID    type_ID;     /* 'adtl' */
} ListHeader;

typedef struct {
	ChunkID    chunk_ID;  /* 'labl' */
	int32_t    chunk_size;

	int32_t    identifier;
	char       text[];
} LabelChunk;


typedef struct {
	ChunkID    chunk_ID; /* 'note' */
	int32_t    chunk_size;

	int32_t    identifier;
	char       text[];
} NoteChunk;


typedef struct {
	ChunkID    chunk_ID; /* 'ltxt' */
	int32_t    chunk_size;

	int32_t    iodentifier;
	int32_t    sample_length;
	int32_t    purpose;
	int16_t    country;
	int16_t    language;
	int16_t    dialect;
	int16_t    code_page;
	char       text[];
} LabelTextChunk;

/* Sampler Chunk */
typedef struct {
	int32_t    identifier;
	int32_t    type;
	int32_t    start;
	int32_t    end;
	int32_t    fraction;
	int32_t    play_count;
} SampleLoop;


typedef struct {
	ChunkID    chunk_ID; /* 'smpl' */
	int32_t    chunk_size;

	int32_t    manufacturer;
	int32_t    product;
	int32_t    sample_period;
	int32_t    midi_unity_note;
	int32_t    MIDI_pitch_fraction;
	int32_t    SMPTE_format;
	int32_t    SMPTE_offset;
	int32_t    sample_loops;
	int32_t    sampler_data;
	SampleLoop loops[];
} SamplerChunk;

/* Instrument Chunk */

typedef struct {
	ChunkID    chunk_ID; /* 'inst' */
	int32_t    chunk_size;

	uint8_t    unshifted_note;
	int8_t     fine_tune;
	int8_t     gain;
	uint8_t    low_note;
	uint8_t    high_note;
	uint8_t    low_velocity;
	uint8_t    high_velocity;
} InstrumentChunk;


static const unsigned char SUB_FORMAT_GUID_PCM[16] = {
	0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x00, 0x80, 0x00, 0x00, 0xAA, 0x00, 0x38, 0x9B, 0x71
};

static const unsigned char SUB_FORMAT_GUID_FLOAT[16] = {
	0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x00, 0x80, 0x00, 0x00, 0xAA, 0x00, 0x38, 0x9B, 0x71
};

#endif /* INTERNAL_RIFFCHUNK_H */
