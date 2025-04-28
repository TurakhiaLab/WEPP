import { b5 as I, b3 as f, aV as q, l as z, aZ as _, aW as E, ar as P } from "./index-8ee3baa4.js";
import { L as v } from "./long-a46ae85c.js";
import { P as l } from "./binary_parser-421d5da5.js";
import "react";
import "react-dom";
const x = 440477507;
function p(c, e) {
  const n = c.prototype[e], t = `_memo_${e}`;
  c.prototype[e] = function() {
    return t in this || (this[t] = n.call(this)), this[t];
  };
}
const g = ["T", "C", "A", "G"], k = [];
for (let c = 0; c < 256; c++)
  k.push(g[c >> 6 & 3] + g[c >> 4 & 3] + g[c >> 2 & 3] + g[c & 3]);
const O = k.map((c) => c.toLowerCase());
class w {
  /**
   * @param {object} args
   * @param {string} [args.path] filesystem path for the .2bit file to open
   * @param {Filehandle} [args.filehandle] node fs.promises-like filehandle for the .2bit file.
   *  Only needs to support `filehandle.read(buffer, offset, length, position)`
   */
  constructor({ filehandle: e, path: n }) {
    if (e)
      this.filehandle = e;
    else if (n)
      this.filehandle = new I(n);
    else
      throw new Error("must supply path or filehandle");
    this.isBigEndian = void 0;
  }
  async _getParser(e) {
    const n = (await this._getParsers())[e];
    if (!n)
      throw new Error(`parser ${e} not found`);
    return n;
  }
  async _detectEndianness() {
    const e = await this.filehandle.read(f.Buffer.allocUnsafe(8), 0, 8, 0), { buffer: n } = e;
    if (n.readInt32LE(0) === x)
      this.isBigEndian = !1, this.version = n.readInt32LE(4);
    else if (n.readInt32BE(0) === x)
      this.isBigEndian = !0, this.version = n.readInt32BE(4);
    else
      throw new Error("not a 2bit file");
  }
  // memoize
  /**
   * @private
   * detects the file's endianness and instantiates our binary parsers accordingly
   */
  async _getParsers() {
    await this._detectEndianness();
    const e = this.isBigEndian ? "big" : "little", n = this.isBigEndian ? "be" : "le";
    let t = new l().endianess(e).uint8("nameLength").string("name", { length: "nameLength" });
    return this.version === 1 ? t = t.buffer("offsetBytes", {
      length: 8
    }) : t = t.uint32("offset"), {
      header: new l().endianess(e).int32("magic", {
        assert: (i) => i === 440477507
      }).int32("version", {
        /* istanbul ignore next */
        assert: (i) => i === 0 || i === 1
      }).uint32("sequenceCount", {
        /* istanbul ignore next */
        assert: (i) => i >= 0
      }).uint32("reserved"),
      index: new l().endianess(e).uint32("sequenceCount").uint32("reserved").array("index", {
        length: "sequenceCount",
        type: t
      }),
      record1: new l().endianess(e).uint32("dnaSize").uint32("nBlockCount"),
      record2: new l().endianess(e).uint32("nBlockCount").array("nBlockStarts", {
        length: "nBlockCount",
        type: `uint32${n}`
      }).array("nBlockSizes", {
        length: "nBlockCount",
        type: `uint32${n}`
      }).uint32("maskBlockCount"),
      record3: new l().endianess(e).uint32("maskBlockCount").array("maskBlockStarts", {
        length: "maskBlockCount",
        type: `uint32${n}`
      }).array("maskBlockSizes", {
        length: "maskBlockCount",
        type: `uint32${n}`
      }).int32("reserved")
      // .buffer('packedDna', { length: 'dnaSize' }),
    };
  }
  // memoize
  /**
   * @returns {Promise} for object with the file's header information, like
   *  `{ magic: 0x1a412743, version: 0, sequenceCount: 42, reserved: 0 }`
   */
  async getHeader() {
    await this._detectEndianness();
    const { buffer: e } = await this.filehandle.read(f.Buffer.allocUnsafe(16), 0, 16, 0);
    return (await this._getParser("header")).parse(e).result;
  }
  // memoize
  /**
   * @returns {Promise} for object with the file's index of offsets, like `{ seqName: fileOffset, ...}`
   */
  async getIndex() {
    const n = 8 + (await this.getHeader()).sequenceCount * (1 + 256 + (this.version === 1 ? 8 : 4)), { buffer: t } = await this.filehandle.read(f.Buffer.allocUnsafe(n), 0, n, 8), i = (await this._getParser("index")).parse(t).result.index, s = {};
    return this.version === 1 ? i.forEach(({ name: o, offsetBytes: a }) => {
      const u = v.fromBytes(a, !0, !this.isBigEndian);
      if (u.greaterThan(Number.MAX_SAFE_INTEGER))
        throw new Error("integer overflow. File offset greater than 2^53-1 encountered. This library can only handle offsets up to 2^53-1.");
      s[o] = u.toNumber();
    }) : i.forEach(({ name: o, offset: a }) => {
      s[o] = a;
    }), s;
  }
  /**
   * @returns {Promise} for an array of string sequence names that are found in the file
   */
  async getSequenceNames() {
    const e = await this.getIndex();
    return Object.keys(e);
  }
  /**
   * @returns {Promise} for an object listing the lengths of all sequences like `{seqName: length, ...}`
   */
  async getSequenceSizes() {
    const e = await this.getIndex(), n = Object.keys(e), t = Object.values(e).map((s) => this._getSequenceSize(s)), r = await Promise.all(t), i = {};
    for (let s = 0; s < n.length; s += 1)
      i[n[s]] = r[s];
    return i;
  }
  /**
   * @param {string} seqName name of the sequence
   * @returns {Promise} for the sequence's length, or undefined if it is not in the file
   */
  async getSequenceSize(e) {
    const t = (await this.getIndex())[e];
    if (t)
      return this._getSequenceSize(t);
  }
  async _getSequenceSize(e) {
    if (e === void 0 || e < 0)
      throw new Error("invalid offset");
    return (await this._parseItem(e, 8, "record1")).dnaSize;
  }
  async _getSequenceRecord(e) {
    if (e === void 0 || e < 0)
      throw new Error("invalid offset");
    const n = await this._parseItem(e, 8, "record1"), t = n.nBlockCount * 8 + 8, r = await this._parseItem(e + 4, t, "record2"), i = r.maskBlockCount * 8 + 8, s = await this._parseItem(e + 4 + t - 4, i, "record3");
    return {
      dnaSize: n.dnaSize,
      nBlocks: { starts: r.nBlockStarts, sizes: r.nBlockSizes },
      maskBlocks: { starts: s.maskBlockStarts, sizes: s.maskBlockSizes },
      dnaPosition: e + 4 + t - 4 + i
    };
  }
  async _parseItem(e, n, t) {
    const { buffer: r } = await this.filehandle.read(f.Buffer.allocUnsafe(n), 0, n, e);
    return (await this._getParser(t)).parse(r).result;
  }
  /**
   * @param {string} seqName name of the sequence you want
   * @param {number} [regionStart] optional 0-based half-open start of the sequence region to fetch.
   * @param {number} [regionEnd] optional 0-based half-open end of the sequence region to fetch. defaults to end of the sequence
   * @returns {Promise} for a string of sequence bases
   */
  async getSequence(e, n = 0, t) {
    const i = (await this.getIndex())[e];
    if (!i)
      return;
    const s = await this._getSequenceRecord(i);
    if (n < 0)
      throw new TypeError("regionStart cannot be less than 0");
    (t === void 0 || t > s.dnaSize) && (t = s.dnaSize);
    const o = this._getOverlappingBlocks(n, t, s.nBlocks.starts, s.nBlocks.sizes), a = this._getOverlappingBlocks(n, t, s.maskBlocks.starts, s.maskBlocks.sizes), u = f.Buffer.allocUnsafe(Math.ceil((t - n) / 4) + 1), h = Math.floor(n / 4), { buffer: C } = await this.filehandle.read(u, 0, u.length, s.dnaPosition + h);
    let m = "";
    for (let d = n; d < t; d += 1) {
      for (; a.length && a[0].end <= d; )
        a.shift();
      const y = a[0] && a[0].start <= d && a[0].end > d;
      if (o[0] && d >= o[0].start && d < o[0].end) {
        const B = o.shift();
        for (; d < B.end && d < t; d += 1)
          m += y ? "n" : "N";
        d -= 1;
      } else {
        const B = Math.floor(d / 4) - h, b = d % 4, S = C[B];
        m += y ? O[S][b] : k[S][b];
      }
    }
    return m;
  }
  _getOverlappingBlocks(e, n, t, r) {
    let i, s;
    for (let a = 0; a < t.length; a += 1) {
      const u = t[a], h = r[a];
      if (e >= u + h || n <= u) {
        if (i !== void 0) {
          s = a;
          break;
        }
      } else
        i === void 0 && (i = a);
    }
    if (i === void 0)
      return [];
    s === void 0 && (s = t.length);
    const o = new Array(s - i);
    for (let a = i; a < s; a += 1)
      o[a - i] = {
        start: t[a],
        end: t[a] + r[a],
        size: r[a]
      };
    return o;
  }
}
p(w, "_getParsers");
p(w, "getIndex");
p(w, "getHeader");
class A extends q.BaseSequenceAdapter {
  async initChromSizes() {
    const e = z.readConfObject(this.config, "chromSizesLocation");
    if (e.uri !== "/path/to/default.chrom.sizes" && e.uri !== "") {
      const t = await _.openLocation(e, this.pluginManager).readFile("utf8");
      return Object.fromEntries(t == null ? void 0 : t.split(/\n|\r\n|\r/).filter((r) => !!r.trim()).map((r) => {
        const [i, s] = r.split("	");
        return [i, +s];
      }));
    }
  }
  constructor(e, n, t) {
    super(e, n, t), this.chromSizesData = this.initChromSizes(), this.twobit = new w({
      filehandle: _.openLocation(z.readConfObject(e, "twoBitLocation"), this.pluginManager)
    });
  }
  async getRefNames() {
    const e = await this.chromSizesData;
    return e ? Object.keys(e) : this.twobit.getSequenceNames();
  }
  async getRegions() {
    const e = await this.chromSizesData;
    if (e)
      return Object.keys(e).map((t) => ({
        refName: t,
        start: 0,
        end: e[t]
      }));
    const n = await this.twobit.getSequenceSizes();
    return Object.keys(n).map((t) => ({
      refName: t,
      start: 0,
      end: n[t]
    }));
  }
  /**
   * Fetch features for a certain region
   * @param param -
   * @returns Observable of Feature objects in the region
   */
  getFeatures({ refName: e, start: n, end: t }) {
    return E(async (r) => {
      const i = await this.twobit.getSequenceSize(e), s = i !== void 0 ? Math.min(i, t) : t, o = await this.twobit.getSequence(e, n, s);
      o && r.next(new P({
        id: `${e} ${n}-${s}`,
        data: { refName: e, start: n, end: s, seq: o }
      })), r.complete();
    });
  }
  /**
   * called to provide a hint that data tied to a certain region
   * will not be needed for the foreseeable future and can be purged
   * from caches, etc
   */
  freeResources() {
  }
}
export {
  A as default
};
