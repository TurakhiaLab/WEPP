import { b9 as z, b5 as N, b8 as Y, b3 as J, bb as K } from "./index-8ee3baa4.js";
import { Q as W } from "./index-05a6222b.js";
import { L as G } from "./long-a46ae85c.js";
function H(a) {
  if (a.greaterThan(Number.MAX_SAFE_INTEGER) || a.lessThan(Number.MIN_SAFE_INTEGER))
    throw new Error("integer overflow");
  return a.toNumber();
}
class Z extends Error {
}
function v(a) {
  if (a && a.aborted) {
    if (typeof DOMException < "u")
      throw new DOMException("aborted", "AbortError");
    {
      const e = new Z("aborted");
      throw e.code = "ERR_ABORTED", e;
    }
  }
}
function ee(a, e) {
  return e.minv.blockPosition - a.maxv.blockPosition < 65e3 && e.maxv.blockPosition - a.minv.blockPosition < 5e6;
}
function j(a, e) {
  const t = [];
  let n = null;
  return a.length === 0 ? a : (a.sort(function(r, i) {
    const o = r.minv.blockPosition - i.minv.blockPosition;
    return o !== 0 ? o : r.minv.dataPosition - i.minv.dataPosition;
  }), a.forEach((r) => {
    (!e || r.maxv.compareTo(e) > 0) && (n === null ? (t.push(r), n = r) : ee(n, r) ? r.maxv.compareTo(n.maxv) > 0 && (n.maxv = r.maxv) : (t.push(r), n = r));
  }), t);
}
class $ {
  constructor(e, t) {
    this.blockPosition = e, this.dataPosition = t;
  }
  toString() {
    return `${this.blockPosition}:${this.dataPosition}`;
  }
  compareTo(e) {
    return this.blockPosition - e.blockPosition || this.dataPosition - e.dataPosition;
  }
  static min(...e) {
    let t, n = 0;
    for (; !t; n += 1)
      t = e[n];
    for (; n < e.length; n += 1)
      t.compareTo(e[n]) > 0 && (t = e[n]);
    return t;
  }
}
function L(a, e = 0, t = !1) {
  if (t)
    throw new Error("big-endian virtual file offsets not implemented");
  return new $(a[e + 7] * 1099511627776 + a[e + 6] * 4294967296 + a[e + 5] * 16777216 + a[e + 4] * 65536 + a[e + 3] * 256 + a[e + 2], a[e + 1] << 8 | a[e]);
}
class P {
  /**
   * @param {VirtualOffset} minv
   * @param {VirtualOffset} maxv
   * @param {number} bin
   * @param {number} [fetchedSize]
   */
  constructor(e, t, n, r = void 0) {
    this.minv = e, this.maxv = t, this.bin = n, this._fetchedSize = r;
  }
  toUniqueString() {
    return `${this.minv}..${this.maxv} (bin ${this.bin}, fetchedSize ${this.fetchedSize()})`;
  }
  toString() {
    return this.toUniqueString();
  }
  compareTo(e) {
    return this.minv.compareTo(e.minv) || this.maxv.compareTo(e.maxv) || this.bin - e.bin;
  }
  fetchedSize() {
    return this._fetchedSize !== void 0 ? this._fetchedSize : this.maxv.blockPosition + 65536 - this.minv.blockPosition;
  }
}
class Q {
  constructor({ filehandle: e, renameRefSeqs: t = (n) => n }) {
    this.filehandle = e, this.renameRefSeq = t;
  }
  async getMetadata(e = {}) {
    const { indices: t, ...n } = await this.parse(e);
    return n;
  }
  _findFirstData(e, t) {
    return e ? e.compareTo(t) > 0 ? t : e : t;
  }
  async parse(e = {}) {
    return this.parseP || (this.parseP = this._parse(e).catch((t) => {
      throw this.parseP = void 0, t;
    })), this.parseP;
  }
  async hasRefSeq(e, t = {}) {
    return !!((await this.parse(t)).indices[e] || {}).binIndex;
  }
}
const te = 21578324, O = 14;
function ne(a, e) {
  return a += 1, e -= 1, [
    [0, 0],
    [1 + (a >> 26), 1 + (e >> 26)],
    [9 + (a >> 23), 9 + (e >> 23)],
    [73 + (a >> 20), 73 + (e >> 20)],
    [585 + (a >> 17), 585 + (e >> 17)],
    [4681 + (a >> 14), 4681 + (e >> 14)]
  ];
}
class A extends Q {
  async lineCount(e, t = {}) {
    const n = await this.parse(t);
    if (!n)
      return -1;
    const r = n.refNameToId[e];
    if (!n.indices[r])
      return -1;
    const { stats: o } = n.indices[r];
    return o ? o.lineCount : -1;
  }
  // fetch and parse the index
  async _parse(e = {}) {
    const t = await this.filehandle.readFile(e), n = await z(t);
    if (v(e.signal), n.readUInt32LE(0) !== te)
      throw new Error("Not a TBI file");
    const r = n.readInt32LE(4), i = n.readInt32LE(8), o = i & 65536 ? "zero-based-half-open" : "1-based-closed", u = {
      0: "generic",
      1: "SAM",
      2: "VCF"
    }[i & 15];
    if (!u)
      throw new Error(`invalid Tabix preset format flags ${i}`);
    const c = {
      ref: n.readInt32LE(12),
      start: n.readInt32LE(16),
      end: n.readInt32LE(20)
    }, h = n.readInt32LE(24), f = 5, l = ((1 << (f + 1) * 3) - 1) / 7, b = 2 ** (14 + f * 3), x = h ? String.fromCharCode(h) : null, w = n.readInt32LE(28), I = n.readInt32LE(32), { refNameToId: m, refIdToName: p } = this._parseNameBytes(n.slice(36, 36 + I));
    let d = 36 + I, g;
    return {
      indices: new Array(r).fill(0).map(() => {
        const E = n.readInt32LE(d);
        d += 4;
        const k = {};
        let R;
        for (let C = 0; C < E; C += 1) {
          const y = n.readUInt32LE(d);
          if (d += 4, y > l + 1)
            throw new Error("tabix index contains too many bins, please use a CSI index");
          if (y === l + 1) {
            const S = n.readInt32LE(d);
            d += 4, S === 2 && (R = this.parsePseudoBin(n, d)), d += 16 * S;
          } else {
            const S = n.readInt32LE(d);
            d += 4;
            const M = new Array(S);
            for (let B = 0; B < S; B += 1) {
              const V = L(n, d), X = L(n, d + 8);
              d += 16, g = this._findFirstData(g, V), M[B] = new P(V, X, y);
            }
            k[y] = M;
          }
        }
        const F = n.readInt32LE(d);
        d += 4;
        const _ = new Array(F);
        for (let C = 0; C < F; C += 1)
          _[C] = L(n, d), d += 8, g = this._findFirstData(g, _[C]);
        return { binIndex: k, linearIndex: _, stats: R };
      }),
      metaChar: x,
      maxBinNumber: l,
      maxRefLength: b,
      skipLines: w,
      firstDataLine: g,
      columnNumbers: c,
      coordinateType: o,
      format: u,
      refIdToName: p,
      refNameToId: m,
      maxBlockSize: 65536
    };
  }
  parsePseudoBin(e, t) {
    return { lineCount: H(G.fromBytesLE(e.slice(t + 16, t + 24), !0)) };
  }
  _parseNameBytes(e) {
    let t = 0, n = 0;
    const r = [], i = {};
    for (let o = 0; o < e.length; o += 1)
      if (!e[o]) {
        if (n < o) {
          let s = e.toString("utf8", n, o);
          s = this.renameRefSeq(s), r[t] = s, i[s] = t;
        }
        n = o + 1, t += 1;
      }
    return { refNameToId: i, refIdToName: r };
  }
  async blocksForRange(e, t, n, r = {}) {
    t < 0 && (t = 0);
    const i = await this.parse(r);
    if (!i)
      return [];
    const o = i.refNameToId[e], s = i.indices[o];
    if (!s)
      return [];
    (s.linearIndex.length ? s.linearIndex[t >> O >= s.linearIndex.length ? s.linearIndex.length - 1 : t >> O] : new $(0, 0)) || console.warn("querying outside of possible tabix range");
    const c = ne(t, n), h = [];
    for (const [w, I] of c)
      for (let m = w; m <= I; m++)
        if (s.binIndex[m]) {
          const p = s.binIndex[m];
          for (let d = 0; d < p.length; ++d)
            h.push(new P(p[d].minv, p[d].maxv, m));
        }
    const f = s.linearIndex.length;
    let l = null;
    const b = Math.min(t >> 14, f - 1), x = Math.min(n >> 14, f - 1);
    for (let w = b; w <= x; ++w) {
      const I = s.linearIndex[w];
      I && (!l || I.compareTo(l) < 0) && (l = I);
    }
    return j(h, l);
  }
}
const ie = 21582659, re = 38359875;
function oe(a, e) {
  return a * 2 ** e;
}
function q(a, e) {
  return Math.floor(a / 2 ** e);
}
class U extends Q {
  constructor(e) {
    super(e), this.maxBinNumber = 0, this.depth = 0, this.minShift = 0;
  }
  async lineCount(e, t = {}) {
    const n = await this.parse(t);
    if (!n)
      return -1;
    const r = n.refNameToId[e];
    if (!n.indices[r])
      return -1;
    const { stats: o } = n.indices[r];
    return o ? o.lineCount : -1;
  }
  async indexCov() {
    throw new Error("CSI indexes do not support indexcov");
  }
  parseAuxData(e, t) {
    const n = e.readInt32LE(t), r = n & 65536 ? "zero-based-half-open" : "1-based-closed", i = { 0: "generic", 1: "SAM", 2: "VCF" }[n & 15];
    if (!i)
      throw new Error(`invalid Tabix preset format flags ${n}`);
    const o = {
      ref: e.readInt32LE(t + 4),
      start: e.readInt32LE(t + 8),
      end: e.readInt32LE(t + 12)
    }, s = e.readInt32LE(t + 16), u = s ? String.fromCharCode(s) : null, c = e.readInt32LE(t + 20), h = e.readInt32LE(t + 24), { refIdToName: f, refNameToId: l } = this._parseNameBytes(e.slice(t + 28, t + 28 + h));
    return {
      refIdToName: f,
      refNameToId: l,
      skipLines: c,
      metaChar: u,
      columnNumbers: o,
      format: i,
      coordinateType: r
    };
  }
  _parseNameBytes(e) {
    let t = 0, n = 0;
    const r = [], i = {};
    for (let o = 0; o < e.length; o += 1)
      if (!e[o]) {
        if (n < o) {
          let s = e.toString("utf8", n, o);
          s = this.renameRefSeq(s), r[t] = s, i[s] = t;
        }
        n = o + 1, t += 1;
      }
    return { refNameToId: i, refIdToName: r };
  }
  // fetch and parse the index
  async _parse(e = {}) {
    const t = await z(await this.filehandle.readFile(e));
    let n;
    if (t.readUInt32LE(0) === ie)
      n = 1;
    else if (t.readUInt32LE(0) === re)
      n = 2;
    else
      throw new Error("Not a CSI file");
    this.minShift = t.readInt32LE(4), this.depth = t.readInt32LE(8), this.maxBinNumber = ((1 << (this.depth + 1) * 3) - 1) / 7;
    const r = 2 ** (this.minShift + this.depth * 3), i = t.readInt32LE(12), o = i && i >= 30 ? this.parseAuxData(t, 16) : {
      refIdToName: [],
      refNameToId: {},
      metaChar: null,
      columnNumbers: { ref: 0, start: 1, end: 2 },
      coordinateType: "zero-based-half-open",
      format: "generic"
    }, s = t.readInt32LE(16 + i);
    let u, c = 16 + i + 4;
    const h = new Array(s).fill(0).map(() => {
      const f = t.readInt32LE(c);
      c += 4;
      const l = {};
      let b;
      for (let x = 0; x < f; x += 1) {
        const w = t.readUInt32LE(c);
        if (w > this.maxBinNumber)
          b = this.parsePseudoBin(t, c + 4), c += 4 + 8 + 4 + 16 + 16;
        else {
          const I = L(t, c + 4);
          u = this._findFirstData(u, I);
          const m = t.readInt32LE(c + 12);
          c += 16;
          const p = new Array(m);
          for (let d = 0; d < m; d += 1) {
            const g = L(t, c), T = L(t, c + 8);
            c += 16, p[d] = new P(g, T, w);
          }
          l[w] = p;
        }
      }
      return { binIndex: l, stats: b };
    });
    return {
      ...o,
      csi: !0,
      refCount: s,
      maxBlockSize: 65536,
      firstDataLine: u,
      csiVersion: n,
      indices: h,
      depth: this.depth,
      maxBinNumber: this.maxBinNumber,
      maxRefLength: r
    };
  }
  parsePseudoBin(e, t) {
    return { lineCount: H(G.fromBytesLE(e.slice(t + 28, t + 36), !0)) };
  }
  async blocksForRange(e, t, n, r = {}) {
    t < 0 && (t = 0);
    const i = await this.parse(r);
    if (!i)
      return [];
    const o = i.refNameToId[e], s = i.indices[o];
    if (!s)
      return [];
    const u = this.reg2bins(t, n), c = [];
    for (const [h, f] of u)
      for (let l = h; l <= f; l++)
        if (s.binIndex[l]) {
          const b = s.binIndex[l];
          for (let x = 0; x < b.length; ++x)
            c.push(new P(b[x].minv, b[x].maxv, l));
        }
    return j(c, new $(0, 0));
  }
  /**
   * calculate the list of bins that may overlap with region [beg,end) (zero-based half-open)
   */
  reg2bins(e, t) {
    e -= 1, e < 1 && (e = 1), t > 2 ** 50 && (t = 2 ** 34), t -= 1;
    let n = 0, r = 0, i = this.minShift + this.depth * 3;
    const o = [];
    for (; n <= this.depth; i -= 3, r += oe(1, n * 3), n += 1) {
      const s = r + q(e, i), u = r + q(t, i);
      if (u - s + o.length > this.maxBinNumber)
        throw new Error(`query ${e}-${t} is too large for current binning scheme (shift ${this.minShift}, depth ${this.depth}), try a smaller query or a coarser index binning scheme`);
      o.push([s, u]);
    }
    return o;
  }
}
const D = typeof TextDecoder < "u" ? new TextDecoder("utf-8") : void 0;
function se(a) {
  return new Promise((e) => setTimeout(e, a));
}
class de {
  /**
   * @param {object} args
   * @param {string} [args.path]
   * @param {filehandle} [args.filehandle]
   * @param {string} [args.tbiPath]
   * @param {filehandle} [args.tbiFilehandle]
   * @param {string} [args.csiPath]
   * @param {filehandle} [args.csiFilehandle]
   * @param {number} [args.yieldTime] yield to main thread after N milliseconds if reading features is taking a long time to avoid hanging main thread
   * @param {function} [args.renameRefSeqs] optional function with sig `string => string` to transform
   * reference sequence names for the purpose of indexing and querying. note that the data that is returned is
   * not altered, just the names of the reference sequences that are used for querying.
   */
  constructor({ path: e, filehandle: t, tbiPath: n, tbiFilehandle: r, csiPath: i, csiFilehandle: o, yieldTime: s = 500, chunkSizeLimit: u = 5e7, renameRefSeqs: c = (f) => f, chunkCacheSize: h = 5 * 2 ** 20 }) {
    if (t)
      this.filehandle = t;
    else if (e)
      this.filehandle = new N(e);
    else
      throw new TypeError("must provide either filehandle or path");
    if (r)
      this.index = new A({
        filehandle: r,
        renameRefSeqs: c
      });
    else if (o)
      this.index = new U({
        filehandle: o,
        renameRefSeqs: c
      });
    else if (n)
      this.index = new A({
        filehandle: new N(n),
        renameRefSeqs: c
      });
    else if (i)
      this.index = new U({
        filehandle: new N(i),
        renameRefSeqs: c
      });
    else if (e)
      this.index = new A({
        filehandle: new N(`${e}.tbi`),
        renameRefSeqs: c
      });
    else
      throw new TypeError("must provide one of tbiFilehandle, tbiPath, csiFilehandle, or csiPath");
    this.chunkSizeLimit = u, this.renameRefSeq = c, this.yieldTime = s, this.chunkCache = new Y({
      cache: new W({ maxSize: Math.floor(h / 65536) }),
      fill: (f, l) => this.readChunk(f, { signal: l })
    });
  }
  /**
   * @param refName name of the reference sequence
   * @param start start of the region (in 0-based half-open coordinates)
   * @param end end of the region (in 0-based half-open coordinates)
   * @param opts callback called for each line in the region. can also pass a object param containing obj.lineCallback, obj.signal, etc
   * @returns promise that is resolved when the whole read is finished, rejected on error
   */
  async getLines(e, t, n, r) {
    let i, o = {}, s;
    if (typeof r > "u")
      throw new TypeError("line callback must be provided");
    if (typeof r == "function" ? s = r : (o = r, s = r.lineCallback), e === void 0)
      throw new TypeError("must provide a reference sequence name");
    if (!s)
      throw new TypeError("line callback must be provided");
    const u = await this.index.getMetadata(o);
    if (v(i), t || (t = 0), n || (n = u.maxRefLength), !(t <= n))
      throw new TypeError("invalid start and end coordinates. start must be less than or equal to end");
    if (t === n)
      return;
    const c = await this.index.blocksForRange(e, t, n, o);
    v(i);
    for (let f = 0; f < c.length; f += 1) {
      const l = c[f].fetchedSize();
      if (l > this.chunkSizeLimit)
        throw new Error(`Too much data. Chunk size ${l.toLocaleString()} bytes exceeds chunkSizeLimit of ${this.chunkSizeLimit.toLocaleString()}.`);
    }
    let h = Date.now();
    for (let f = 0; f < c.length; f += 1) {
      let l;
      const b = c[f], { buffer: x, cpositions: w, dpositions: I } = await this.chunkCache.get(b.toString(), b);
      v(i);
      let m = 0, p = 0;
      for (; m < x.length; ) {
        const d = x.indexOf(`
`, m);
        if (d === -1)
          break;
        const g = x.slice(m, d), T = (D == null ? void 0 : D.decode(g)) || g.toString();
        if (I) {
          for (; m + b.minv.dataPosition >= I[p++]; )
            ;
          p--;
        }
        const { startCoordinate: E, overlaps: k } = this.checkLine(u, e, t, n, T);
        if (l !== void 0 && E !== void 0 && l > E)
          throw new Error(`Lines not sorted by start coordinate (${l} > ${E}), this file is not usable with Tabix.`);
        if (l = E, k)
          s(
            T.trim(),
            // cpositions[pos] refers to actual file offset of a bgzip block boundaries
            //
            // we multiply by (1 <<8) in order to make sure each block has a "unique"
            // address space so that data in that block could never overlap
            //
            // then the blockStart-dpositions is an uncompressed file offset from
            // that bgzip block boundary, and since the cpositions are multiplied by
            // (1 << 8) these uncompressed offsets get a unique space
            w[p] * 256 + (m - I[p]) + b.minv.dataPosition + 1
          );
        else if (E !== void 0 && E >= n)
          return;
        this.yieldTime && h - Date.now() > this.yieldTime && (h = Date.now(), v(i), await se(1)), m = d + 1;
      }
    }
  }
  async getMetadata(e = {}) {
    return this.index.getMetadata(e);
  }
  /**
   * get a buffer containing the "header" region of
   * the file, which are the bytes up to the first
   * non-meta line
   */
  async getHeaderBuffer(e = {}) {
    const { firstDataLine: t, metaChar: n, maxBlockSize: r } = await this.getMetadata(e);
    v(e.signal);
    const i = ((t == null ? void 0 : t.blockPosition) || 0) + r;
    let o = await this._readRegion(0, i, e);
    v(e.signal);
    try {
      o = await z(o);
    } catch (s) {
      throw console.error(s), new Error(
        //@ts-ignore
        `error decompressing block ${s.code} at 0 (length ${i}) ${s}`
      );
    }
    if (n) {
      let s = -1;
      const u = `
`.charCodeAt(0), c = n.charCodeAt(0);
      for (let h = 0; h < o.length && !(h === s + 1 && o[h] !== c); h += 1)
        o[h] === u && (s = h);
      o = o.slice(0, s + 1);
    }
    return o;
  }
  /**
   * get a string containing the "header" region of the
   * file, is the portion up to the first non-meta line
   *
   * @returns {Promise} for a string
   */
  async getHeader(e = {}) {
    return (await this.getHeaderBuffer(e)).toString("utf8");
  }
  /**
   * get an array of reference sequence names, in the order in which
   * they occur in the file. reference sequence renaming is not applied
   * to these names.
   */
  async getReferenceSequenceNames(e = {}) {
    return (await this.getMetadata(e)).refIdToName;
  }
  /**
   * @param {object} metadata metadata object from the parsed index,
   * containing columnNumbers, metaChar, and format
   * @param {string} regionRefName
   * @param {number} regionStart region start coordinate (0-based-half-open)
   * @param {number} regionEnd region end coordinate (0-based-half-open)
   * @param {array[string]} line
   * @returns {object} like `{startCoordinate, overlaps}`. overlaps is boolean,
   * true if line is a data line that overlaps the given region
   */
  checkLine(e, t, n, r, i) {
    const { columnNumbers: o, metaChar: s, coordinateType: u, format: c } = e;
    if (i.charAt(0) === s)
      return { overlaps: !1 };
    let { ref: h, start: f, end: l } = o;
    h || (h = 0), f || (f = 0), l || (l = 0), c === "VCF" && (l = 8);
    const b = Math.max(h, f, l);
    let x = 1, w = 0, I = "", m = -1 / 0;
    for (let p = 0; p < i.length + 1; p += 1)
      if (i[p] === "	" || p === i.length) {
        if (x === h) {
          if (this.renameRefSeq(i.slice(w, p)) !== t)
            return { overlaps: !1 };
        } else if (x === f) {
          if (m = parseInt(i.slice(w, p), 10), u === "1-based-closed" && (m -= 1), m >= r)
            return { startCoordinate: m, overlaps: !1 };
          if ((l === 0 || l === f) && m + 1 <= n)
            return { startCoordinate: m, overlaps: !1 };
        } else if (c === "VCF" && x === 4)
          I = i.slice(w, p);
        else if (x === l) {
          let d;
          if (c === "VCF" ? d = this._getVcfEnd(m, I, i.slice(w, p)) : d = parseInt(i.slice(w, p), 10), d <= n)
            return { overlaps: !1 };
        }
        if (w = p + 1, x += 1, x > b)
          break;
      }
    return { startCoordinate: m, overlaps: !0 };
  }
  _getVcfEnd(e, t, n) {
    let r = e + t.length;
    const i = n.indexOf("SVTYPE=TRA") !== -1;
    if (n[0] !== "." && !i) {
      let o = ";";
      for (let s = 0; s < n.length; s += 1) {
        if (o === ";" && n.slice(s, s + 4) === "END=") {
          let u = n.indexOf(";", s);
          u === -1 && (u = n.length), r = parseInt(n.slice(s + 4, u), 10);
          break;
        }
        o = n[s];
      }
    } else if (i)
      return e + 1;
    return r;
  }
  /**
   * return the approximate number of data lines in the given reference sequence
   * @param refSeq reference sequence name
   * @returns number of data lines present on that reference sequence
   */
  async lineCount(e, t = {}) {
    return this.index.lineCount(e, t);
  }
  async _readRegion(e, t, n = {}) {
    const r = J.Buffer.alloc(t), { bytesRead: i, buffer: o } = await this.filehandle.read(r, 0, t, e, n);
    return o.slice(0, i);
  }
  /**
   * read and uncompress the data in a chunk (composed of one or more
   * contiguous bgzip blocks) of the file
   */
  async readChunk(e, t = {}) {
    const n = await this._readRegion(e.minv.blockPosition, e.fetchedSize(), t);
    try {
      return K(n, e);
    } catch (r) {
      throw new Error(`error decompressing c ${e.toString()} ${r}`);
    }
  }
}
export {
  de as T
};
