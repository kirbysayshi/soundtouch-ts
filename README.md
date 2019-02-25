## SoundTouch-TS

A port of a port to TypeScript. Pitchshift and Timestretch in JS/TS. This library is LGPL2.1 due to [SoundTouch](https://gitlab.com/soundtouch/soundtouch). A more tested port of this library, with more utilities, is available in [soundtouch-js](https://github.com/cutterbl/SoundTouchJS). This TS port exists because I wanted the types, and didn't know soundtouch-js existed until recently.

However, as long as you allow a user, at runtime, to swap out this library, it should fall under section 6b of https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html. But JS/TS and LGPL are legally ambiguous, so use at your own risk.

## Usage

```ts
// It's TS, but JS is fine too!
import { SoundTouch } from 'soundtouch-ts';
const st = new SoundTouch(44100);

const tempo = 0.5;

// Audio will take 2x as long to play with no pitch changes
st.tempo = tempo

const ctx = new AudioContext;

fetch('http://test-audio.somewhere.mp3')
.then(res => res.arrayBuffer())
.then(ab => ctx.decodeAudio(ab))
.then(ab => {
  const interleaved = asInterleaved(ab);
  st.inputBuffer.putSamples(interleaved);
  st.process();

  const channels = 2;
  const receiver = new Float32Array(channels * ab.length / tempo)

  const requested = 256;
  let received = 0;
  while (st.outputBuffer.frameCount) {
    const queued = st.frameCount;
    st.process();
    st.outputBuffer.receiveSamples(receiver.subarray(received), requested);
    const remaining = st.outputBuffer.frameCount;
    received = (queued - remaining) * channels;
  }

  const audio = asPlanar(receiver, 44100, 2);
  const abnode = new AudioBufferSourceNode(ctx, { buffer: audio });
  ctx.destination.connect(abnode);
  abnode.start();
})

function asInterleaved(ab: AudioBuffer): Float32Array {
  const channels = ab.numberOfChannels;
  const output = new Float32Array(channels * ab.length);
  for (let i = 0; i < ab.length; i++) {
    for (let c = 0; c < channels; c++) {
      const chan = ab.getChannelData(c);
      output[i * channels + c] = chan[i];
    }
  }
  return output;
}

function asPlanar(
  buffer: Float32Array,
  sampleRate: number,
  channels: number = 2
): AudioBuffer {
  const channelLength = Math.floor(buffer.length / channels);
  const output = new AudioBuffer({
    numberOfChannels: channels,
    length: channelLength,
    sampleRate: sampleRate
  });

  for (let c = 0; c < channels; c++) {
    const chan = output.getChannelData(c);
    for (let i = 0; i < channelLength; i++) {
      chan[i] = buffer[i * channels + c];
    }
  }

  return output;
}
```

## Publishing

```sh
$ npx pack publish
```

## History

This port was modified from the following:

- Original Port (LGPL 2.1): https://github.com/also/soundtouch-js/tree/master/src/js
- Modularized / Expanded (MIT): https://github.com/jakubfiala/soundtouch-js
- Modified and included in a UI (MIT): https://github.com/ZVK/stretcher (http://zackzukowski.com/TAP-audio-player/)
- Converted to TypeScript: [src/index.ts](src/index.ts).

## License

[GNU Lesser General Public Library, version 2.1](https://www.gnu.org/licenses/lgpl-2.1.en.html)
