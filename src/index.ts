/*
 * SoundTouch JS audio processing library
 * Copyright (c) Olli Parviainen
 * Copyright (c) Ryan Berdeen
 * Copyright (c) Zach Zukowski
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/**
 * Giving this value for the sequence length sets automatic parameter value
 * according to tempo setting (recommended)
 */
const USE_AUTO_SEQUENCE_LEN = 0;

/**
 * Default length of a single processing sequence, in milliseconds. This determines to how
 * long sequences the original sound is chopped in the time-stretch algorithm.
 *
 * The larger this value is, the lesser sequences are used in processing. In principle
 * a bigger value sounds better when slowing down tempo, but worse when increasing tempo
 * and vice versa.
 *
 * Increasing this value reduces computational burden and vice versa.
 */
//var DEFAULT_SEQUENCE_MS = 130
const DEFAULT_SEQUENCE_MS = USE_AUTO_SEQUENCE_LEN;

/**
 * Giving this value for the seek window length sets automatic parameter value
 * according to tempo setting (recommended)
 */
const USE_AUTO_SEEKWINDOW_LEN = 0;

/**
 * Seeking window default length in milliseconds for algorithm that finds the best possible
 * overlapping location. This determines from how wide window the algorithm may look for an
 * optimal joining location when mixing the sound sequences back together.
 *
 * The bigger this window setting is, the higher the possibility to find a better mixing
 * position will become, but at the same time large values may cause a "drifting" artifact
 * because consequent sequences will be taken at more uneven intervals.
 *
 * If there's a disturbing artifact that sounds as if a constant frequency was drifting
 * around, try reducing this setting.
 *
 * Increasing this value increases computational burden and vice versa.
 */
//var DEFAULT_SEEKWINDOW_MS = 25;
const DEFAULT_SEEKWINDOW_MS = USE_AUTO_SEEKWINDOW_LEN;

/**
 * Overlap length in milliseconds. When the chopped sound sequences are mixed back together,
 * to form a continuous sound stream, this parameter defines over how long period the two
 * consecutive sequences are let to overlap each other.
 *
 * This shouldn't be that critical parameter. If you reduce the DEFAULT_SEQUENCE_MS setting
 * by a large amount, you might wish to try a smaller value on this.
 *
 * Increasing this value increases computational burden and vice versa.
 */
const DEFAULT_OVERLAP_MS = 8;

// Table for the hierarchical mixing position seeking algorithm
// prettier-ignore
const _SCAN_OFFSETS = [
  [ 124,  186,  248,  310,  372,  434,  496,  558,  620,  682,  744, 806,
    868,  930,  992, 1054, 1116, 1178, 1240, 1302, 1364, 1426, 1488,   0],
  [-100,  -75,  -50,  -25,   25,   50,   75,  100,    0,    0,    0,   0,
      0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,   0],
  [ -20,  -15,  -10,   -5,    5,   10,   15,   20,    0,    0,    0,   0,
      0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,   0],
  [  -4,   -3,   -2,   -1,    1,    2,    3,    4,    0,    0,    0,   0,
      0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,   0]];

// Adjust tempo param according to tempo, so that variating processing sequence length is used
// at varius tempo settings, between the given low...top limits
const AUTOSEQ_TEMPO_LOW = 0.5; // auto setting low tempo range (-50%)
const AUTOSEQ_TEMPO_TOP = 2.0; // auto setting top tempo range (+100%)

// sequence-ms setting values at above low & top tempo
const AUTOSEQ_AT_MIN = 125.0;
const AUTOSEQ_AT_MAX = 50.0;
const AUTOSEQ_K =
  (AUTOSEQ_AT_MAX - AUTOSEQ_AT_MIN) / (AUTOSEQ_TEMPO_TOP - AUTOSEQ_TEMPO_LOW);
const AUTOSEQ_C = AUTOSEQ_AT_MIN - AUTOSEQ_K * AUTOSEQ_TEMPO_LOW;

// seek-window-ms setting values at above low & top tempo
const AUTOSEEK_AT_MIN = 25.0;
const AUTOSEEK_AT_MAX = 15.0;
const AUTOSEEK_K =
  (AUTOSEEK_AT_MAX - AUTOSEEK_AT_MIN) / (AUTOSEQ_TEMPO_TOP - AUTOSEQ_TEMPO_LOW);
const AUTOSEEK_C = AUTOSEEK_AT_MIN - AUTOSEEK_K * AUTOSEQ_TEMPO_LOW;

function testFloatEqual(a: number, b: number) {
  return (a > b ? a - b : b - a) > 1e-10;
}

class AbstractFifoSamplePipe {
  inputBuffer: FifoSampleBuffer | null;
  outputBuffer: FifoSampleBuffer | null;

  constructor(createBuffers: boolean) {
    if (createBuffers) {
      this.inputBuffer = new FifoSampleBuffer();
      this.outputBuffer = new FifoSampleBuffer();
    } else {
      this.inputBuffer = this.outputBuffer = null;
    }
  }

  clear() {
    if (this.inputBuffer) this.inputBuffer.clear();
    if (this.outputBuffer) this.outputBuffer.clear();
  }
}

class RateTransposer extends AbstractFifoSamplePipe {
  rate: number;

  slopeCount: number = 0;
  prevSampleL: number = 0;
  prevSampleR: number = 0;

  constructor(createBuffers: boolean) {
    super(createBuffers);
    this._reset();
    this.rate = 1;
  }

  _reset() {
    this.slopeCount = 0;
    this.prevSampleL = 0;
    this.prevSampleR = 0;
  }

  process() {
    // TODO aa filter
    const numFrames = this.inputBuffer!.frameCount;
    this.outputBuffer!.ensureAdditionalCapacity(numFrames / this.rate + 1);
    const numFramesOutput = this._transpose(numFrames);
    this.inputBuffer!.receive();
    this.outputBuffer!.put(numFramesOutput);
  }

  _transpose(numFrames: number) {
    if (numFrames === 0) {
      return 0; // No work.
    }

    const src = this.inputBuffer!.vector;
    const srcOffset = this.inputBuffer!.startIndex;

    const dest = this.outputBuffer!.vector;
    const destOffset = this.outputBuffer!.endIndex;

    let used = 0;
    let i = 0;

    while (this.slopeCount < 1.0) {
      dest[destOffset + 2 * i] =
        (1.0 - this.slopeCount) * this.prevSampleL +
        this.slopeCount * src[srcOffset];
      dest[destOffset + 2 * i + 1] =
        (1.0 - this.slopeCount) * this.prevSampleR +
        this.slopeCount * src[srcOffset + 1];
      i++;
      this.slopeCount += this.rate;
    }

    this.slopeCount -= 1.0;

    if (numFrames != 1) {
      out: while (true) {
        while (this.slopeCount > 1.0) {
          this.slopeCount -= 1.0;
          used++;
          if (used >= numFrames - 1) {
            break out;
          }
        }

        const srcIndex = srcOffset + 2 * used;
        dest[destOffset + 2 * i] =
          (1.0 - this.slopeCount) * src[srcIndex] +
          this.slopeCount * src[srcIndex + 2];
        dest[destOffset + 2 * i + 1] =
          (1.0 - this.slopeCount) * src[srcIndex + 1] +
          this.slopeCount * src[srcIndex + 3];

        i++;
        this.slopeCount += this.rate;
      }
    }

    this.prevSampleL = src[srcOffset + 2 * numFrames - 2];
    this.prevSampleR = src[srcOffset + 2 * numFrames - 1];

    return i;
  }
}

class FifoSampleBuffer {
  private _vector = new Float32Array(0);
  private _position = 0;
  private _frameCount = 0;

  constructor() {}

  get vector() {
    return this._vector;
  }

  get position() {
    return this._position;
  }

  get startIndex() {
    return this._position * 2;
  }

  get frameCount() {
    return this._frameCount;
  }

  get endIndex() {
    return (this._position + this._frameCount) * 2;
  }

  clear(frameCount: number = -1) {
    this.receive(frameCount);
    this.rewind();
  }

  put(numFrames: number) {
    this._frameCount += numFrames;
  }

  putSamples(samples: Float32Array, position = 0, numFrames = -1) {
    const sourceOffset = position * 2;
    if (!(numFrames >= 0)) {
      numFrames = (samples.length - sourceOffset) / 2;
    }
    const numSamples = numFrames * 2;

    this.ensureCapacity(numFrames + this._frameCount);

    const destOffset = this.endIndex;
    this._vector.set(
      samples.subarray(sourceOffset, sourceOffset + numSamples),
      destOffset
    );

    this._frameCount += numFrames;
  }

  putBuffer(buffer: FifoSampleBuffer, position = 0, numFrames = -1) {
    if (!(numFrames >= 0)) {
      numFrames = buffer.frameCount - position;
    }
    this.putSamples(buffer.vector, buffer.position + position, numFrames);
  }

  receive(numFrames: number = -1) {
    if (!(numFrames >= 0) || numFrames > this._frameCount) {
      numFrames = this._frameCount;
    }
    this._frameCount -= numFrames;
    this._position += numFrames;
  }

  receiveSamples(output: Float32Array, numFrames: number) {
    const numSamples = numFrames * 2;
    const sourceOffset = this.startIndex;
    output.set(this._vector.subarray(sourceOffset, sourceOffset + numSamples));
    this.receive(numFrames);
  }

  extract(output: Float32Array, position: number, numFrames: number) {
    const sourceOffset = this.startIndex + position * 2;
    const numSamples = numFrames * 2;
    output.set(this._vector.subarray(sourceOffset, sourceOffset + numSamples));
  }

  ensureCapacity(numFrames: number) {
    const minLength = numFrames * 2;
    if (this._vector.length < minLength) {
      const newVector = new Float32Array(minLength);
      newVector.set(this._vector.subarray(this.startIndex, this.endIndex));
      this._vector = newVector;
      this._position = 0;
    } else {
      this.rewind();
    }
  }

  ensureAdditionalCapacity(numFrames: number) {
    this.ensureCapacity(this.frameCount + numFrames);
  }

  rewind() {
    if (this._position > 0) {
      this._vector.set(this._vector.subarray(this.startIndex, this.endIndex));
      this._position = 0;
    }
  }
}

class Stretch extends AbstractFifoSamplePipe {
  sampleRate: number = 44100;

  sequenceMs: number = DEFAULT_SEQUENCE_MS;
  seekWindowMs: number = DEFAULT_SEEKWINDOW_MS;
  overlapMs: number = DEFAULT_OVERLAP_MS;

  bQuickSeek: boolean = true;
  bMidBufferDirty: boolean = false;

  pRefMidBuffer = new Float32Array(0);
  pMidBuffer: Float32Array | null = new Float32Array(0);
  overlapLength: number = 0;

  bAutoSeqSetting: boolean = true;
  bAutoSeekSetting: boolean = true;

  nominalSkip: number = 0;
  skipFract: number = 0;

  seekWindowLength: number = 0;
  seekLength: number = 0;
  sampleReq: number = 0;

  private _tempo: number = 1;

  constructor(createBuffers: boolean, sampleRate: number) {
    super(createBuffers);

    this.setParameters(
      sampleRate,
      DEFAULT_SEQUENCE_MS,
      DEFAULT_SEEKWINDOW_MS,
      DEFAULT_OVERLAP_MS
    );
  }

  clear() {
    AbstractFifoSamplePipe.prototype.clear.call(this);
    this._clearMidBuffer();
  }

  _clearMidBuffer() {
    if (this.bMidBufferDirty) {
      this.bMidBufferDirty = false;
      this.pMidBuffer = null;
    }
  }

  /**
   * Sets routine control parameters. These control are certain time constants
   * defining how the sound is stretched to the desired duration.
   *
   * 'sampleRate' = sample rate of the sound
   * 'sequenceMS' = one processing sequence length in milliseconds (default = 82 ms)
   * 'seekwindowMS' = seeking window length for scanning the best overlapping
   *      position (default = 28 ms)
   * 'overlapMS' = overlapping length (default = 12 ms)
   */
  setParameters(
    aSampleRate: number,
    aSequenceMS: number,
    aSeekWindowMS: number,
    aOverlapMS: number
  ) {
    // accept only positive parameter values - if zero or negative, use old values instead
    if (aSampleRate > 0) {
      this.sampleRate = aSampleRate;
    }
    if (aOverlapMS > 0) {
      this.overlapMs = aOverlapMS;
    }

    if (aSequenceMS > 0) {
      this.sequenceMs = aSequenceMS;
      this.bAutoSeqSetting = false;
    } else {
      // zero or below, use automatic setting
      this.bAutoSeqSetting = true;
    }

    if (aSeekWindowMS > 0) {
      this.seekWindowMs = aSeekWindowMS;
      this.bAutoSeekSetting = false;
    } else {
      // zero or below, use automatic setting
      this.bAutoSeekSetting = true;
    }

    this.calcSeqParameters();

    this.calculateOverlapLength(this.overlapMs);

    // set tempo to recalculate 'sampleReq'
    this.tempo = this._tempo;
  }

  /**
   * Sets new target tempo. Normal tempo = 'SCALE', smaller values represent slower
   * tempo, larger faster tempo.
   */
  set tempo(newTempo) {
    let intskip;

    this._tempo = newTempo;

    // Calculate new sequence duration
    this.calcSeqParameters();

    // Calculate ideal skip length (according to tempo value)
    this.nominalSkip =
      this._tempo * (this.seekWindowLength - this.overlapLength);
    this.skipFract = 0;
    intskip = Math.floor(this.nominalSkip + 0.5);

    // Calculate how many samples are needed in the 'inputBuffer' to
    // process another batch of samples
    this.sampleReq =
      Math.max(intskip + this.overlapLength, this.seekWindowLength) +
      this.seekLength;
  }

  get tempo() {
    return this._tempo;
  }

  get inputChunkSize() {
    return this.sampleReq;
  }

  get outputChunkSize() {
    return (
      this.overlapLength +
      Math.max(0, this.seekWindowLength - 2 * this.overlapLength)
    );
  }

  /**
   * Calculates overlapInMsec period length in samples.
   */
  calculateOverlapLength(overlapInMsec: number) {
    let newOvl;

    // TODO assert(overlapInMsec >= 0);
    newOvl = (this.sampleRate * overlapInMsec) / 1000;
    if (newOvl < 16) newOvl = 16;

    // must be divisible by 8
    newOvl -= newOvl % 8;

    this.overlapLength = newOvl;

    this.pRefMidBuffer = new Float32Array(this.overlapLength * 2);
    this.pMidBuffer = new Float32Array(this.overlapLength * 2);
  }

  checkLimits(x: number, mi: number, ma: number) {
    return x < mi ? mi : x > ma ? ma : x;
  }

  /**
   * Calculates processing sequence length according to tempo setting
   */
  calcSeqParameters() {
    let seq;
    let seek;

    if (this.bAutoSeqSetting) {
      seq = AUTOSEQ_C + AUTOSEQ_K * this._tempo;
      seq = this.checkLimits(seq, AUTOSEQ_AT_MAX, AUTOSEQ_AT_MIN);
      this.sequenceMs = Math.floor(seq + 0.5);
    }

    if (this.bAutoSeekSetting) {
      seek = AUTOSEEK_C + AUTOSEEK_K * this._tempo;
      seek = this.checkLimits(seek, AUTOSEEK_AT_MAX, AUTOSEEK_AT_MIN);
      this.seekWindowMs = Math.floor(seek + 0.5);
    }

    // Update seek window lengths
    this.seekWindowLength = Math.floor(
      (this.sampleRate * this.sequenceMs) / 1000
    );
    this.seekLength = Math.floor((this.sampleRate * this.seekWindowMs) / 1000);
  }

  /**
   * Enables/disables the quick position seeking algorithm.
   */
  set quickSeek(enable: boolean) {
    this.bQuickSeek = enable;
  }

  /**
   * Seeks for the optimal overlap-mixing position.
   */
  seekBestOverlapPosition() {
    if (this.bQuickSeek) {
      return this.seekBestOverlapPositionStereoQuick();
    } else {
      return this.seekBestOverlapPositionStereo();
    }
  }

  /**
   * Seeks for the optimal overlap-mixing position. The 'stereo' version of the
   * routine
   *
   * The best position is determined as the position where the two overlapped
   * sample sequences are 'most alike', in terms of the highest cross-correlation
   * value over the overlapping period
   */
  seekBestOverlapPositionStereo() {
    let bestOffs;
    let bestCorr;
    let corr;
    let i;

    // Slopes the amplitudes of the 'midBuffer' samples.
    this.precalcCorrReferenceStereo();

    bestCorr = Number.MIN_VALUE;
    bestOffs = 0;

    // Scans for the best correlation value by testing each possible position
    // over the permitted range.
    for (i = 0; i < this.seekLength; i++) {
      // Calculates correlation value for the mixing position corresponding
      // to 'i'
      corr = this.calcCrossCorrStereo(2 * i, this.pRefMidBuffer);

      // Checks for the highest correlation value.
      if (corr > bestCorr) {
        bestCorr = corr;
        bestOffs = i;
      }
    }
    return bestOffs;
  }

  /**
   * Seeks for the optimal overlap-mixing position. The 'stereo' version of the
   * routine
   *
   * The best position is determined as the position where the two overlapped
   * sample sequences are 'most alike', in terms of the highest cross-correlation
   * value over the overlapping period
   */
  seekBestOverlapPositionStereoQuick() {
    let j;
    let bestOffs;
    let bestCorr;
    let corr;
    let scanCount;
    let corrOffset;
    let tempOffset;

    // Slopes the amplitude of the 'midBuffer' samples
    this.precalcCorrReferenceStereo();

    bestCorr = Number.MIN_VALUE;
    bestOffs = 0;
    corrOffset = 0;
    tempOffset = 0;

    // Scans for the best correlation value using four-pass hierarchical search.
    //
    // The look-up table 'scans' has hierarchical position adjusting steps.
    // In first pass the routine searhes for the highest correlation with
    // relatively coarse steps, then rescans the neighbourhood of the highest
    // correlation with better resolution and so on.
    for (scanCount = 0; scanCount < 4; scanCount++) {
      j = 0;
      while (_SCAN_OFFSETS[scanCount][j]) {
        tempOffset = corrOffset + _SCAN_OFFSETS[scanCount][j];
        if (tempOffset >= this.seekLength) {
          break;
        }

        // Calculates correlation value for the mixing position corresponding
        // to 'tempOffset'
        corr = this.calcCrossCorrStereo(2 * tempOffset, this.pRefMidBuffer);

        // Checks for the highest correlation value
        if (corr > bestCorr) {
          bestCorr = corr;
          bestOffs = tempOffset;
        }
        j++;
      }
      corrOffset = bestOffs;
    }
    return bestOffs;
  }

  /**
   * Slopes the amplitude of the 'midBuffer' samples so that cross correlation
   * is faster to calculate
   */
  precalcCorrReferenceStereo() {
    let i;
    let cnt2;
    let temp;

    for (i = 0; i < this.overlapLength; i++) {
      temp = i * (this.overlapLength - i);
      cnt2 = i * 2;
      this.pRefMidBuffer[cnt2] = this.pMidBuffer![cnt2] * temp;
      this.pRefMidBuffer[cnt2 + 1] = this.pMidBuffer![cnt2 + 1] * temp;
    }
  }

  calcCrossCorrStereo(mixingPos: number, compare: Float32Array) {
    const mixing = this.inputBuffer!.vector;
    mixingPos += this.inputBuffer!.startIndex;

    let corr;
    let i;
    let mixingOffset;
    corr = 0;
    for (i = 2; i < 2 * this.overlapLength; i += 2) {
      mixingOffset = i + mixingPos;
      corr +=
        mixing[mixingOffset] * compare[i] +
        mixing[mixingOffset + 1] * compare[i + 1];
    }
    return corr;
  }

  // TODO inline
  /**
   * Overlaps samples in 'midBuffer' with the samples in 'pInputBuffer' at position
   * of 'ovlPos'.
   */
  overlap(ovlPos: number) {
    this.overlapStereo(2 * ovlPos);
  }

  /**
   * Overlaps samples in 'midBuffer' with the samples in 'pInput'
   */
  overlapStereo(pInputPos: number) {
    const pInput = this.inputBuffer!.vector;
    pInputPos += this.inputBuffer!.startIndex;

    const pOutput = this.outputBuffer!.vector;
    const pOutputPos = this.outputBuffer!.endIndex;
    let i;
    let cnt2;
    let fTemp;
    let fScale;
    let fi;
    let pInputOffset;
    let pOutputOffset;

    fScale = 1 / this.overlapLength;
    for (i = 0; i < this.overlapLength; i++) {
      fTemp = (this.overlapLength - i) * fScale;
      fi = i * fScale;
      cnt2 = 2 * i;
      pInputOffset = cnt2 + pInputPos;
      pOutputOffset = cnt2 + pOutputPos;
      pOutput[pOutputOffset + 0] =
        pInput[pInputOffset + 0] * fi + this.pMidBuffer![cnt2 + 0] * fTemp;
      pOutput[pOutputOffset + 1] =
        pInput[pInputOffset + 1] * fi + this.pMidBuffer![cnt2 + 1] * fTemp;
    }
  }

  process() {
    let ovlSkip;
    let offset;
    let temp;
    let i;
    if (this.pMidBuffer === null) {
      // if midBuffer is empty, move the first samples of the input stream
      // into it
      if (this.inputBuffer!.frameCount < this.overlapLength) {
        // wait until we've got overlapLength samples
        return;
      }
      this.pMidBuffer = new Float32Array(this.overlapLength * 2);
      this.inputBuffer!.receiveSamples(this.pMidBuffer, this.overlapLength);
    }

    let output;
    // Process samples as long as there are enough samples in 'inputBuffer'
    // to form a processing frame.
    while (this.inputBuffer!.frameCount >= this.sampleReq) {
      // If tempo differs from the normal ('SCALE'), scan for the best overlapping
      // position
      offset = this.seekBestOverlapPosition();

      // Mix the samples in the 'inputBuffer' at position of 'offset' with the
      // samples in 'midBuffer' using sliding overlapping
      // ... first partially overlap with the end of the previous sequence
      // (that's in 'midBuffer')
      this.outputBuffer!.ensureAdditionalCapacity(this.overlapLength);
      // FIXME unit?
      //overlap(uint(offset));
      this.overlap(Math.floor(offset));
      this.outputBuffer!.put(this.overlapLength);

      // ... then copy sequence samples from 'inputBuffer' to output
      temp = this.seekWindowLength - 2 * this.overlapLength; // & 0xfffffffe;
      if (temp > 0) {
        this.outputBuffer!.putBuffer(
          this.inputBuffer!,
          offset + this.overlapLength,
          temp
        );
      }

      // Copies the end of the current sequence from 'inputBuffer' to
      // 'midBuffer' for being mixed with the beginning of the next
      // processing sequence and so on
      //assert(offset + seekWindowLength <= (int)inputBuffer.numSamples());
      const start =
        this.inputBuffer!.startIndex +
        2 * (offset + this.seekWindowLength - this.overlapLength);
      this.pMidBuffer.set(
        this.inputBuffer!.vector.subarray(start, start + 2 * this.overlapLength)
      );

      // Remove the processed samples from the input buffer. Update
      // the difference between integer & nominal skip step to 'skipFract'
      // in order to prevent the error from accumulating over time.
      this.skipFract += this.nominalSkip; // real skip size
      ovlSkip = Math.floor(this.skipFract); // rounded to integer skip
      this.skipFract -= ovlSkip; // maintain the fraction part, i.e. real vs. integer skip
      this.inputBuffer!.receive(ovlSkip);
    }
  }
}

class SoundTouch {
  rateTransposer: RateTransposer;
  tdStretch: Stretch;

  _inputBuffer: FifoSampleBuffer;
  _intermediateBuffer: FifoSampleBuffer;
  _outputBuffer: FifoSampleBuffer;

  _rate = 0;
  _tempo = 0;

  virtualPitch = 1.0;
  virtualRate = 1.0;
  virtualTempo = 1.0;

  constructor(sampleRate: number) {
    this.rateTransposer = new RateTransposer(false);
    this.tdStretch = new Stretch(false, sampleRate);

    this._inputBuffer = new FifoSampleBuffer();
    this._intermediateBuffer = new FifoSampleBuffer();
    this._outputBuffer = new FifoSampleBuffer();

    this._rate = 0;
    this._tempo = 0;

    this.virtualPitch = 1.0;
    this.virtualRate = 1.0;
    this.virtualTempo = 1.0;

    this._calculateEffectiveRateAndTempo();
  }

  clear() {
    this.rateTransposer.clear();
    this.tdStretch.clear();
  }

  get rate() {
    return this._rate;
  }

  set rate(rate) {
    this.virtualRate = rate;
    this._calculateEffectiveRateAndTempo();
  }

  rateChange(rateChange: number) {
    this.rate = 1.0 + 0.01 * rateChange;
  }

  get tempo() {
    return this._tempo;
  }

  set tempo(tempo) {
    this.virtualTempo = tempo;
    this._calculateEffectiveRateAndTempo();
  }

  set tempoChange(tempoChange: number) {
    this.tempo = 1.0 + 0.01 * tempoChange;
  }

  set pitch(pitch: number) {
    this.virtualPitch = pitch;
    this._calculateEffectiveRateAndTempo();
  }

  set pitchOctaves(pitchOctaves: number) {
    this.pitch = Math.exp(0.69314718056 * pitchOctaves);
    this._calculateEffectiveRateAndTempo();
  }

  set pitchSemitones(pitchSemitones: number) {
    this.pitchOctaves = pitchSemitones / 12.0;
  }

  get inputBuffer() {
    return this._inputBuffer!;
  }

  get outputBuffer() {
    return this._outputBuffer!;
  }

  _calculateEffectiveRateAndTempo() {
    const previousTempo = this._tempo;
    const previousRate = this._rate;

    this._tempo = this.virtualTempo / this.virtualPitch;
    this._rate = this.virtualRate * this.virtualPitch;

    if (testFloatEqual(this._tempo, previousTempo)) {
      this.tdStretch.tempo = this._tempo;
    }
    if (testFloatEqual(this._rate, previousRate)) {
      this.rateTransposer.rate = this._rate;
    }

    if (this._rate > 1.0) {
      if (this.outputBuffer! != this.rateTransposer.outputBuffer) {
        this.tdStretch.inputBuffer = this.inputBuffer!;
        this.tdStretch.outputBuffer = this._intermediateBuffer;

        this.rateTransposer.inputBuffer = this._intermediateBuffer;
        this.rateTransposer.outputBuffer = this.outputBuffer!;
      }
    } else {
      if (this.outputBuffer! != this.tdStretch.outputBuffer) {
        this.rateTransposer.inputBuffer = this.inputBuffer!;
        this.rateTransposer.outputBuffer = this._intermediateBuffer;

        this.tdStretch.inputBuffer = this._intermediateBuffer;
        this.tdStretch.outputBuffer = this.outputBuffer!;
      }
    }
  }

  process() {
    if (this._rate > 1.0) {
      this.tdStretch.process();
      this.rateTransposer.process();
    } else {
      this.rateTransposer.process();
      this.tdStretch.process();
    }
  }
}

export { RateTransposer, Stretch, SoundTouch };
