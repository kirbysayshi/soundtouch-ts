import { SoundTouch } from ".";

test("Basic API", () => {
  const st = new SoundTouch(44100);
  expect(st).toBeTruthy();
  expect(st.virtualPitch).toBe(1);
});
