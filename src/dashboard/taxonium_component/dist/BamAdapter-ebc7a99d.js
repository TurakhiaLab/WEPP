import { L as Y } from "./long-a46ae85c.js";
import { b9 as H, Q as ne, b8 as se, b5 as M, b6 as D, b3 as O, ba as ae, bb as oe, b7 as ce, a1 as he, bc as fe, bd as de, aV as ue, aZ as N, t as V, aX as le, aY as me, aW as ge } from "./index-8ee3baa4.js";
class J {
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
    let t, i = 0;
    for (; !t; i += 1)
      t = e[i];
    for (; i < e.length; i += 1)
      t.compareTo(e[i]) > 0 && (t = e[i]);
    return t;
  }
}
function R(l, e = 0, t = !1) {
  if (t)
    throw new Error("big-endian virtual file offsets not implemented");
  return new J(l[e + 7] * 1099511627776 + l[e + 6] * 4294967296 + l[e + 5] * 16777216 + l[e + 4] * 65536 + l[e + 3] * 256 + l[e + 2], l[e + 1] << 8 | l[e]);
}
class k {
  /**
   * @param {VirtualOffset} minv
   * @param {VirtualOffset} maxv
   * @param {number} bin
   * @param {number} [fetchedSize]
   */
  constructor(e, t, i, r = void 0) {
    this.minv = e, this.maxv = t, this.bin = i, this._fetchedSize = r;
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
class W {
  /**
   * @param {filehandle} filehandle
   * @param {function} [renameRefSeqs]
   */
  constructor({ filehandle: e, renameRefSeq: t = (i) => i }) {
    this.filehandle = e, this.renameRefSeq = t;
  }
  _findFirstData(e, t) {
    const i = e.firstDataLine;
    i ? e.firstDataLine = i.compareTo(t) > 0 ? t : i : e.firstDataLine = t;
  }
  async parse(e = {}) {
    return this.setupP || (this.setupP = this._parse(e).catch((t) => {
      throw this.setupP = void 0, t;
    })), this.setupP;
  }
  async hasRefSeq(e, t = {}) {
    return !!((await this.parse(t)).indices[e] || {}).binIndex;
  }
}
function pe(l) {
  return new Promise((e) => setTimeout(e, l));
}
function K(l) {
  if (l.greaterThan(Number.MAX_SAFE_INTEGER) || l.lessThan(Number.MIN_SAFE_INTEGER))
    throw new Error("integer overflow");
  return l.toNumber();
}
function Z(l) {
  if (l && l.aborted) {
    if (typeof DOMException < "u")
      throw new DOMException("aborted", "AbortError");
    {
      const e = new Error("aborted");
      throw e.code = "ERR_ABORTED", e;
    }
  }
}
async function ee(l) {
  await Promise.resolve(), Z(l);
}
function _e(l, e) {
  return e.minv.blockPosition - l.maxv.blockPosition < 65e3 && e.maxv.blockPosition - l.minv.blockPosition < 5e6;
}
function ye(l = {}) {
  return "aborted" in l ? { signal: l } : l;
}
function te(l, e) {
  const t = [];
  let i = null;
  return l.length === 0 ? l : (l.sort((r, a) => {
    const o = r.minv.blockPosition - a.minv.blockPosition;
    return o !== 0 ? o : r.minv.dataPosition - a.minv.dataPosition;
  }), l.forEach((r) => {
    (!e || r.maxv.compareTo(e) > 0) && (i === null ? (t.push(r), i = r) : _e(i, r) ? r.maxv.compareTo(i.maxv) > 0 && (i.maxv = r.maxv) : (t.push(r), i = r));
  }), t);
}
const be = 21578050;
function we(l, e) {
  return l - l % e;
}
function Ee(l, e) {
  return l - l % e + e;
}
class T extends W {
  parsePseudoBin(e, t) {
    return { lineCount: K(Y.fromBytesLE(Array.prototype.slice.call(e, t + 16, t + 24), !0)) };
  }
  async lineCount(e, t = {}) {
    const r = (await this.parse(t)).indices[e];
    if (!r)
      return -1;
    const a = r.stats || {};
    return a.lineCount === void 0 ? -1 : a.lineCount;
  }
  fetchBai(e = {}) {
    return this.baiP || (this.baiP = this.filehandle.readFile(e).catch((t) => {
      throw this.baiP = void 0, t;
    })), this.baiP;
  }
  // fetch and parse the index
  async _parse() {
    const e = { bai: !0, maxBlockSize: 65536 }, t = await this.fetchBai();
    if (t.readUInt32LE(0) !== be)
      throw new Error("Not a BAI file");
    e.refCount = t.readInt32LE(4);
    const r = ((1 << (5 + 1) * 3) - 1) / 7;
    e.indices = new Array(e.refCount);
    let a = 8;
    for (let o = 0; o < e.refCount; o += 1) {
      const h = t.readInt32LE(a);
      let d;
      a += 4;
      const s = {};
      for (let u = 0; u < h; u += 1) {
        const c = t.readUInt32LE(a);
        if (a += 4, c === r + 1)
          a += 4, d = this.parsePseudoBin(t, a), a += 32;
        else {
          if (c > r + 1)
            throw new Error("bai index contains too many bins, please use CSI");
          {
            const g = t.readInt32LE(a);
            a += 4;
            const _ = new Array(g);
            for (let y = 0; y < g; y += 1) {
              const w = R(t, a), A = R(t, a + 8);
              a += 16, this._findFirstData(e, w), _[y] = new k(w, A, c);
            }
            s[c] = _;
          }
        }
      }
      const m = t.readInt32LE(a);
      a += 4;
      const p = new Array(m);
      for (let u = 0; u < m; u += 1)
        p[u] = R(t, a), a += 8, this._findFirstData(e, p[u]);
      e.indices[o] = { binIndex: s, linearIndex: p, stats: d };
    }
    return e;
  }
  async indexCov(e, t, i, r = {}) {
    const o = t !== void 0, d = (await this.parse(r)).indices[e];
    if (!d)
      return [];
    const { linearIndex: s = [], stats: m } = d;
    if (!s.length)
      return [];
    const p = i !== void 0 ? Ee(i, 16384) : (s.length - 1) * 16384, u = t !== void 0 ? we(t, 16384) : 0;
    let c;
    o ? c = new Array((p - u) / 16384) : c = new Array(s.length - 1);
    const g = s[s.length - 1].blockPosition;
    if (p > (s.length - 1) * 16384)
      throw new Error("query outside of range of linear index");
    let _ = s[u / 16384].blockPosition;
    for (let y = u / 16384, w = 0; y < p / 16384; y++, w++)
      c[w] = {
        score: s[y + 1].blockPosition - _,
        start: y * 16384,
        end: y * 16384 + 16384
      }, _ = s[y + 1].blockPosition;
    return c.map((y) => ({ ...y, score: y.score * m.lineCount / g }));
  }
  /**
   * calculate the list of bins that may overlap with region [beg,end) (zero-based half-open)
   * @returns {Array[number]}
   */
  reg2bins(e, t) {
    return t -= 1, [
      [0, 0],
      [1 + (e >> 26), 1 + (t >> 26)],
      [9 + (e >> 23), 9 + (t >> 23)],
      [73 + (e >> 20), 73 + (t >> 20)],
      [585 + (e >> 17), 585 + (t >> 17)],
      [4681 + (e >> 14), 4681 + (t >> 14)]
    ];
  }
  async blocksForRange(e, t, i, r = {}) {
    t < 0 && (t = 0);
    const a = await this.parse(r);
    if (!a)
      return [];
    const o = a.indices[e];
    if (!o)
      return [];
    const h = this.reg2bins(t, i), d = [];
    for (const [c, g] of h)
      for (let _ = c; _ <= g; _++)
        if (o.binIndex[_]) {
          const y = o.binIndex[_];
          for (let w = 0; w < y.length; ++w)
            d.push(new k(y[w].minv, y[w].maxv, _));
        }
    const s = o.linearIndex.length;
    let m = null;
    const p = Math.min(t >> 14, s - 1), u = Math.min(i >> 14, s - 1);
    for (let c = p; c <= u; ++c) {
      const g = o.linearIndex[c];
      g && (!m || g.compareTo(m) < 0) && (m = g);
    }
    return te(d, m);
  }
}
const xe = 21582659, Ae = 38359875;
function Ie(l, e) {
  return l * 2 ** e;
}
function Q(l, e) {
  return Math.floor(l / 2 ** e);
}
class U extends W {
  constructor(e) {
    super(e), this.maxBinNumber = 0, this.depth = 0, this.minShift = 0;
  }
  async lineCount(e) {
    const t = await this.parse();
    if (!t || !t.indices[e])
      return -1;
    const { stats: r } = t.indices[e];
    return r ? r.lineCount : -1;
  }
  async indexCov() {
    return [];
  }
  parseAuxData(e, t, i) {
    if (i < 30)
      return {};
    const r = {};
    if (r.formatFlags = e.readInt32LE(t), r.coordinateType = r.formatFlags & 65536 ? "zero-based-half-open" : "1-based-closed", r.format = { 0: "generic", 1: "SAM", 2: "VCF" }[r.formatFlags & 15], !r.format)
      throw new Error(`invalid Tabix preset format flags ${r.formatFlags}`);
    r.columnNumbers = {
      ref: e.readInt32LE(t + 4),
      start: e.readInt32LE(t + 8),
      end: e.readInt32LE(t + 12)
    }, r.metaValue = e.readInt32LE(t + 16), r.metaChar = r.metaValue ? String.fromCharCode(r.metaValue) : "", r.skipLines = e.readInt32LE(t + 20);
    const a = e.readInt32LE(t + 24);
    return Object.assign(r, this._parseNameBytes(e.subarray(t + 28, t + 28 + a))), r;
  }
  _parseNameBytes(e) {
    let t = 0, i = 0;
    const r = [], a = {};
    for (let o = 0; o < e.length; o += 1)
      if (!e[o]) {
        if (i < o) {
          let h = e.toString("utf8", i, o);
          h = this.renameRefSeq(h), r[t] = h, a[h] = t;
        }
        i = o + 1, t += 1;
      }
    return { refNameToId: a, refIdToName: r };
  }
  // fetch and parse the index
  async _parse(e) {
    const t = { csi: !0, maxBlockSize: 65536 }, i = await this.filehandle.readFile(e), r = await H(i);
    if (r.readUInt32LE(0) === xe)
      t.csiVersion = 1;
    else if (r.readUInt32LE(0) === Ae)
      t.csiVersion = 2;
    else
      throw new Error("Not a CSI file");
    this.minShift = r.readInt32LE(4), this.depth = r.readInt32LE(8), this.maxBinNumber = ((1 << (this.depth + 1) * 3) - 1) / 7;
    const a = r.readInt32LE(12);
    a && Object.assign(t, this.parseAuxData(r, 16, a)), t.refCount = r.readInt32LE(16 + a), t.indices = new Array(t.refCount);
    let o = 16 + a + 4;
    for (let h = 0; h < t.refCount; h += 1) {
      await ee(e.signal);
      const d = r.readInt32LE(o);
      o += 4;
      const s = {};
      let m;
      for (let p = 0; p < d; p += 1) {
        const u = r.readUInt32LE(o);
        if (u > this.maxBinNumber)
          m = this.parsePseudoBin(r, o + 4), o += 4 + 8 + 4 + 16 + 16;
        else {
          const c = R(r, o + 4);
          this._findFirstData(t, c);
          const g = r.readInt32LE(o + 12);
          o += 16;
          const _ = new Array(g);
          for (let y = 0; y < g; y += 1) {
            const w = R(r, o), A = R(r, o + 8);
            o += 16, _[y] = new k(w, A, u);
          }
          s[u] = _;
        }
      }
      t.indices[h] = { binIndex: s, stats: m };
    }
    return t;
  }
  parsePseudoBin(e, t) {
    return { lineCount: K(Y.fromBytesLE(Array.prototype.slice.call(e, t + 28, t + 36), !0)) };
  }
  async blocksForRange(e, t, i, r = {}) {
    t < 0 && (t = 0);
    const a = await this.parse(r), o = a == null ? void 0 : a.indices[e];
    if (!o)
      return [];
    const h = this.reg2bins(t, i), d = [];
    for (const [s, m] of h)
      for (let p = s; p <= m; p++)
        if (o.binIndex[p]) {
          const u = o.binIndex[p];
          for (let c = 0; c < u.length; ++c)
            d.push(new k(u[c].minv, u[c].maxv, p));
        }
    return te(d, new J(0, 0));
  }
  /**
   * calculate the list of bins that may overlap with region [beg,end)
   * (zero-based half-open)
   */
  reg2bins(e, t) {
    e -= 1, e < 1 && (e = 1), t > 2 ** 50 && (t = 2 ** 34), t -= 1;
    let i = 0, r = 0, a = this.minShift + this.depth * 3;
    const o = [];
    for (; i <= this.depth; a -= 3, r += Ie(1, i * 3), i += 1) {
      const h = r + Q(e, a), d = r + Q(t, a);
      if (d - h + o.length > this.maxBinNumber)
        throw new Error(`query ${e}-${t} is too large for current binning scheme (shift ${this.minShift}, depth ${this.depth}), try a smaller query or a coarser index binning scheme`);
      o.push([h, d]);
    }
    return o;
  }
}
class Se {
  constructor(e = {}) {
    if (!(e.maxSize && e.maxSize > 0))
      throw new TypeError("`maxSize` must be a number greater than 0");
    this.maxSize = e.maxSize, this.cache = /* @__PURE__ */ new Map(), this.oldCache = /* @__PURE__ */ new Map(), this._size = 0;
  }
  _set(e, t) {
    this.cache.set(e, t), this._size++, this._size >= this.maxSize && (this._size = 0, this.oldCache = this.cache, this.cache = /* @__PURE__ */ new Map());
  }
  get(e) {
    if (this.cache.has(e))
      return this.cache.get(e);
    if (this.oldCache.has(e)) {
      const t = this.oldCache.get(e);
      return this._set(e, t), t;
    }
  }
  set(e, t) {
    return this.cache.has(e) ? this.cache.set(e, t) : this._set(e, t), this;
  }
  has(e) {
    return this.cache.has(e) || this.oldCache.has(e);
  }
  peek(e) {
    if (this.cache.has(e))
      return this.cache.get(e);
    if (this.oldCache.has(e))
      return this.oldCache.get(e);
  }
  delete(e) {
    const t = this.cache.delete(e);
    return t && this._size--, this.oldCache.delete(e) || t;
  }
  clear() {
    this.cache.clear(), this.oldCache.clear(), this._size = 0;
  }
  *keys() {
    for (const [e] of this)
      yield e;
  }
  *values() {
    for (const [, e] of this)
      yield e;
  }
  *[Symbol.iterator]() {
    for (const e of this.cache)
      yield e;
    for (const e of this.oldCache) {
      const [t] = e;
      this.cache.has(t) || (yield e);
    }
  }
  get size() {
    let e = 0;
    for (const t of this.oldCache.keys())
      this.cache.has(t) || e++;
    return this._size + e;
  }
}
var Ce = Se;
const ve = /* @__PURE__ */ ne(Ce), C = {
  //  the read is paired in sequencing, no matter whether it is mapped in a pair
  BAM_FPAIRED: 1,
  //  the read is mapped in a proper pair
  BAM_FPROPER_PAIR: 2,
  //  the read itself is unmapped; conflictive with BAM_FPROPER_PAIR
  BAM_FUNMAP: 4,
  //  the mate is unmapped
  BAM_FMUNMAP: 8,
  //  the read is mapped to the reverse strand
  BAM_FREVERSE: 16,
  //  the mate is mapped to the reverse strand
  BAM_FMREVERSE: 32,
  //  this is read1
  BAM_FREAD1: 64,
  //  this is read2
  BAM_FREAD2: 128,
  //  not primary alignment
  BAM_FSECONDARY: 256,
  //  QC failure
  BAM_FQCFAIL: 512,
  //  optical or PCR duplicate
  BAM_FDUP: 1024,
  //  supplementary alignment
  BAM_FSUPPLEMENTARY: 2048
}, X = "=ACMGRSVTWYHKDBN".split(""), B = "MIDNSHP=X???????".split("");
class Le {
  constructor(e) {
    this.data = {}, this._tagList = [], this._allTagsParsed = !1;
    const { bytes: t, fileOffset: i } = e, { byteArray: r, start: a } = t;
    this.data = {}, this.bytes = t, this._id = i, this._refID = r.readInt32LE(a + 4), this.data.start = r.readInt32LE(a + 8), this.flags = (r.readInt32LE(a + 16) & 4294901760) >> 16;
  }
  get(e) {
    return this[e] ? this.data[e] ? this.data[e] : (this.data[e] = this[e](), this.data[e]) : this._get(e.toLowerCase());
  }
  end() {
    return this.get("start") + this.get("length_on_ref");
  }
  seq_id() {
    return this._refID;
  }
  // same as get(), except requires lower-case arguments.  used
  // internally to save lots of calls to field.toLowerCase()
  _get(e) {
    return e in this.data ? this.data[e] : (this.data[e] = this._parseTag(e), this.data[e]);
  }
  _tags() {
    this._parseAllTags();
    let e = ["seq"];
    this.isSegmentUnmapped() || e.push("start", "end", "strand", "score", "qual", "MQ", "CIGAR", "length_on_ref", "template_length"), this.isPaired() && e.push("next_segment_position", "pair_orientation"), e = e.concat(this._tagList || []), Object.keys(this.data).forEach((i) => {
      i[0] !== "_" && i !== "next_seq_id" && e.push(i);
    });
    const t = {};
    return e.filter((i) => {
      if (i in this.data && this.data[i] === void 0 || i === "CG" || i === "cg")
        return !1;
      const r = i.toLowerCase(), a = t[r];
      return t[r] = !0, !a;
    });
  }
  parent() {
  }
  children() {
    return this.get("subfeatures");
  }
  id() {
    return this._id;
  }
  // special parsers
  /**
   * Mapping quality score.
   */
  mq() {
    const e = (this.get("_bin_mq_nl") & 65280) >> 8;
    return e === 255 ? void 0 : e;
  }
  score() {
    return this.get("mq");
  }
  qual() {
    var e;
    return (e = this.qualRaw()) === null || e === void 0 ? void 0 : e.join(" ");
  }
  qualRaw() {
    if (this.isSegmentUnmapped())
      return;
    const { start: e, byteArray: t } = this.bytes, i = e + 36 + this.get("_l_read_name") + this.get("_n_cigar_op") * 4 + this.get("_seq_bytes"), r = this.get("seq_length");
    return t.subarray(i, i + r);
  }
  strand() {
    return this.isReverseComplemented() ? -1 : 1;
  }
  multi_segment_next_segment_strand() {
    if (!this.isMateUnmapped())
      return this.isMateReverseComplemented() ? -1 : 1;
  }
  name() {
    return this.get("_read_name");
  }
  _read_name() {
    const e = this.get("_l_read_name"), { byteArray: t, start: i } = this.bytes;
    return t.toString("ascii", i + 36, i + 36 + e - 1);
  }
  /**
   * Get the value of a tag, parsing the tags as far as necessary.
   * Only called if we have not already parsed that field.
   */
  _parseTag(e) {
    if (this._allTagsParsed)
      return;
    const { byteArray: t, start: i } = this.bytes;
    let r = this._tagOffset || i + 36 + this.get("_l_read_name") + this.get("_n_cigar_op") * 4 + this.get("_seq_bytes") + this.get("seq_length");
    const a = this.bytes.end;
    let o;
    for (; r < a && o !== e; ) {
      const h = String.fromCharCode(t[r], t[r + 1]);
      o = h.toLowerCase();
      const d = String.fromCharCode(t[r + 2]);
      r += 3;
      let s;
      switch (d) {
        case "A":
          s = String.fromCharCode(t[r]), r += 1;
          break;
        case "i":
          s = t.readInt32LE(r), r += 4;
          break;
        case "I":
          s = t.readUInt32LE(r), r += 4;
          break;
        case "c":
          s = t.readInt8(r), r += 1;
          break;
        case "C":
          s = t.readUInt8(r), r += 1;
          break;
        case "s":
          s = t.readInt16LE(r), r += 2;
          break;
        case "S":
          s = t.readUInt16LE(r), r += 2;
          break;
        case "f":
          s = t.readFloatLE(r), r += 4;
          break;
        case "Z":
        case "H":
          for (s = ""; r <= a; ) {
            const m = t[r++];
            if (m === 0)
              break;
            s += String.fromCharCode(m);
          }
          break;
        case "B": {
          s = "";
          const m = t[r++], p = String.fromCharCode(m), u = t.readInt32LE(r);
          if (r += 4, p === "i")
            if (h === "CG")
              for (let c = 0; c < u; c++) {
                const g = t.readInt32LE(r), _ = g >> 4, y = B[g & 15];
                s += _ + y, r += 4;
              }
            else
              for (let c = 0; c < u; c++)
                s += t.readInt32LE(r), c + 1 < u && (s += ","), r += 4;
          if (p === "I")
            if (h === "CG")
              for (let c = 0; c < u; c++) {
                const g = t.readUInt32LE(r), _ = g >> 4, y = B[g & 15];
                s += _ + y, r += 4;
              }
            else
              for (let c = 0; c < u; c++)
                s += t.readUInt32LE(r), c + 1 < u && (s += ","), r += 4;
          if (p === "s")
            for (let c = 0; c < u; c++)
              s += t.readInt16LE(r), c + 1 < u && (s += ","), r += 2;
          if (p === "S")
            for (let c = 0; c < u; c++)
              s += t.readUInt16LE(r), c + 1 < u && (s += ","), r += 2;
          if (p === "c")
            for (let c = 0; c < u; c++)
              s += t.readInt8(r), c + 1 < u && (s += ","), r += 1;
          if (p === "C")
            for (let c = 0; c < u; c++)
              s += t.readUInt8(r), c + 1 < u && (s += ","), r += 1;
          if (p === "f")
            for (let c = 0; c < u; c++)
              s += t.readFloatLE(r), c + 1 < u && (s += ","), r += 4;
          break;
        }
        default:
          console.warn(`Unknown BAM tag type '${d}', tags may be incomplete`), s = void 0, r = a;
      }
      if (this._tagOffset = r, this._tagList.push(h), o === e)
        return s;
      this.data[o] = s;
    }
    this._allTagsParsed = !0;
  }
  _parseAllTags() {
    this._parseTag("");
  }
  _parseCigar(e) {
    return (
      //@ts-ignore
      e.match(/\d+\D/g).map((t) => [t.match(/\D/)[0].toUpperCase(), parseInt(t, 10)])
    );
  }
  /**
   * @returns {boolean} true if the read is paired, regardless of whether both segments are mapped
   */
  isPaired() {
    return !!(this.flags & C.BAM_FPAIRED);
  }
  /** @returns {boolean} true if the read is paired, and both segments are mapped */
  isProperlyPaired() {
    return !!(this.flags & C.BAM_FPROPER_PAIR);
  }
  /** @returns {boolean} true if the read itself is unmapped; conflictive with isProperlyPaired */
  isSegmentUnmapped() {
    return !!(this.flags & C.BAM_FUNMAP);
  }
  /** @returns {boolean} true if the read itself is unmapped; conflictive with isProperlyPaired */
  isMateUnmapped() {
    return !!(this.flags & C.BAM_FMUNMAP);
  }
  /** @returns {boolean} true if the read is mapped to the reverse strand */
  isReverseComplemented() {
    return !!(this.flags & C.BAM_FREVERSE);
  }
  /** @returns {boolean} true if the mate is mapped to the reverse strand */
  isMateReverseComplemented() {
    return !!(this.flags & C.BAM_FMREVERSE);
  }
  /** @returns {boolean} true if this is read number 1 in a pair */
  isRead1() {
    return !!(this.flags & C.BAM_FREAD1);
  }
  /** @returns {boolean} true if this is read number 2 in a pair */
  isRead2() {
    return !!(this.flags & C.BAM_FREAD2);
  }
  /** @returns {boolean} true if this is a secondary alignment */
  isSecondary() {
    return !!(this.flags & C.BAM_FSECONDARY);
  }
  /** @returns {boolean} true if this read has failed QC checks */
  isFailedQc() {
    return !!(this.flags & C.BAM_FQCFAIL);
  }
  /** @returns {boolean} true if the read is an optical or PCR duplicate */
  isDuplicate() {
    return !!(this.flags & C.BAM_FDUP);
  }
  /** @returns {boolean} true if this is a supplementary alignment */
  isSupplementary() {
    return !!(this.flags & C.BAM_FSUPPLEMENTARY);
  }
  cigar() {
    if (this.isSegmentUnmapped())
      return;
    const { byteArray: e, start: t } = this.bytes, i = this.get("_n_cigar_op");
    let r = t + 36 + this.get("_l_read_name");
    const a = this.get("seq_length");
    let o = "", h = 0, d = e.readInt32LE(r), s = d >> 4, m = B[d & 15];
    if (m === "S" && s === a)
      return r += 4, d = e.readInt32LE(r), s = d >> 4, m = B[d & 15], m !== "N" && console.warn("CG tag with no N tag"), this.data.length_on_ref = s, this.get("CG");
    for (let p = 0; p < i; ++p)
      d = e.readInt32LE(r), s = d >> 4, m = B[d & 15], o += s + m, m !== "H" && m !== "S" && m !== "I" && (h += s), r += 4;
    return this.data.length_on_ref = h, o;
  }
  _flags() {
  }
  length_on_ref() {
    return this.data.length_on_ref ? this.data.length_on_ref : (this.get("cigar"), this.data.length_on_ref);
  }
  _n_cigar_op() {
    return this.get("_flag_nc") & 65535;
  }
  _l_read_name() {
    return this.get("_bin_mq_nl") & 255;
  }
  /**
   * number of bytes in the sequence field
   */
  _seq_bytes() {
    return this.get("seq_length") + 1 >> 1;
  }
  getReadBases() {
    return this.seq();
  }
  seq() {
    const { byteArray: e, start: t } = this.bytes, i = t + 36 + this.get("_l_read_name") + this.get("_n_cigar_op") * 4, r = this.get("_seq_bytes"), a = this.get("seq_length");
    let o = "", h = 0;
    for (let d = 0; d < r; ++d) {
      const s = e[i + d];
      o += X[(s & 240) >> 4], h++, h < a && (o += X[s & 15], h++);
    }
    return o;
  }
  // adapted from igv.js
  getPairOrientation() {
    if (!this.isSegmentUnmapped() && !this.isMateUnmapped() && this._refID === this._next_refid()) {
      const e = this.isReverseComplemented() ? "R" : "F", t = this.isMateReverseComplemented() ? "R" : "F";
      let i = " ", r = " ";
      this.isRead1() ? (i = "1", r = "2") : this.isRead2() && (i = "2", r = "1");
      const a = [];
      return this.template_length() > 0 ? (a[0] = e, a[1] = i, a[2] = t, a[3] = r) : (a[2] = e, a[3] = i, a[0] = t, a[1] = r), a.join("");
    }
    return null;
  }
  _bin_mq_nl() {
    return this.bytes.byteArray.readInt32LE(this.bytes.start + 12);
  }
  _flag_nc() {
    return this.bytes.byteArray.readInt32LE(this.bytes.start + 16);
  }
  seq_length() {
    return this.bytes.byteArray.readInt32LE(this.bytes.start + 20);
  }
  _next_refid() {
    return this.bytes.byteArray.readInt32LE(this.bytes.start + 24);
  }
  _next_pos() {
    return this.bytes.byteArray.readInt32LE(this.bytes.start + 28);
  }
  template_length() {
    return this.bytes.byteArray.readInt32LE(this.bytes.start + 32);
  }
  toJSON() {
    const e = {};
    return Object.keys(this).forEach((t) => {
      t.charAt(0) === "_" || t === "bytes" || (e[t] = this[t]);
    }), e;
  }
}
function Re(l) {
  const e = l.split(/\r?\n/), t = [];
  return e.forEach((i) => {
    const [r, ...a] = i.split(/\t/), o = a.map((h) => {
      const [d, s] = h.split(":", 2);
      return { tag: d, value: s };
    });
    r && t.push({ tag: r.substr(1), data: o });
  }), t;
}
const Pe = 21840194, $ = 65536;
function j(l) {
  return [].concat(...l);
}
async function Te(l) {
  const e = [];
  for await (const t of l)
    e.push(t);
  return e;
}
class Be {
  /**
   * @param {object} args
   * @param {string} [args.bamPath]
   * @param {FileHandle} [args.bamFilehandle]
   * @param {string} [args.baiPath]
   * @param {FileHandle} [args.baiFilehandle]
   */
  constructor({ bamFilehandle: e, bamPath: t, bamUrl: i, baiPath: r, baiFilehandle: a, baiUrl: o, csiPath: h, csiFilehandle: d, csiUrl: s, fetchSizeLimit: m, chunkSizeLimit: p, yieldThreadTime: u = 100, renameRefSeqs: c = (g) => g }) {
    if (this.featureCache = new se({
      //@ts-ignore
      cache: new ve({
        maxSize: 50
      }),
      //@ts-ignore
      fill: async ({ chunk: g, opts: _ }, y) => {
        const { data: w, cpositions: A, dpositions: I } = await this._readChunk({
          chunk: g,
          opts: { ..._, signal: y }
        });
        return await this.readBamFeatures(w, A, I, g);
      }
    }), this.renameRefSeq = c, e)
      this.bam = e;
    else if (t)
      this.bam = new M(t);
    else if (i)
      this.bam = new D(i);
    else
      throw new Error("unable to initialize bam");
    if (d)
      this.index = new U({ filehandle: d });
    else if (h)
      this.index = new U({ filehandle: new M(h) });
    else if (s)
      this.index = new U({ filehandle: new D(s) });
    else if (a)
      this.index = new T({ filehandle: a });
    else if (r)
      this.index = new T({ filehandle: new M(r) });
    else if (o)
      this.index = new T({ filehandle: new D(o) });
    else if (t)
      this.index = new T({ filehandle: new M(`${t}.bai`) });
    else if (i)
      this.index = new T({ filehandle: new D(`${i}.bai`) });
    else
      throw new Error("unable to infer index format");
    this.fetchSizeLimit = m || 5e8, this.chunkSizeLimit = p || 3e8, this.yieldThreadTime = u;
  }
  async getHeader(e = {}) {
    const t = ye(e), i = await this.index.parse(t), r = i.firstDataLine ? i.firstDataLine.blockPosition + 65535 : void 0;
    let a;
    if (r) {
      const m = await this.bam.read(O.Buffer.alloc(r + $), 0, r + $, 0, t), { bytesRead: p } = m;
      if ({ buffer: a } = m, !p)
        throw new Error("Error reading header");
      p < r ? a = a.subarray(0, p) : a = a.subarray(0, r);
    } else
      a = await this.bam.readFile(t);
    const o = await H(a);
    if (o.readInt32LE(0) !== Pe)
      throw new Error("Not a BAM file");
    const h = o.readInt32LE(4);
    this.header = o.toString("utf8", 8, 8 + h);
    const { chrToIndex: d, indexToChr: s } = await this._readRefSeqs(h + 8, 65535, t);
    return this.chrToIndex = d, this.indexToChr = s, Re(this.header);
  }
  async getHeaderText(e = {}) {
    return await this.getHeader(e), this.header;
  }
  // the full length of the refseq block is not given in advance so this grabs
  // a chunk and doubles it if all refseqs haven't been processed
  async _readRefSeqs(e, t, i = {}) {
    if (e > t)
      return this._readRefSeqs(e, t * 2, i);
    const r = t + $, { bytesRead: a, buffer: o } = await this.bam.read(O.Buffer.alloc(r), 0, t, 0, i);
    if (!a)
      throw new Error("Error reading refseqs from header");
    const h = await H(o.subarray(0, Math.min(a, t))), d = h.readInt32LE(e);
    let s = e + 4;
    const m = {}, p = [];
    for (let u = 0; u < d; u += 1) {
      const c = h.readInt32LE(s), g = this.renameRefSeq(h.toString("utf8", s + 4, s + 4 + c - 1)), _ = h.readInt32LE(s + c + 4);
      if (m[g] = u, p.push({ refName: g, length: _ }), s = s + 8 + c, s > h.length)
        return console.warn(`BAM header is very big.  Re-fetching ${t} bytes.`), this._readRefSeqs(e, t * 2, i);
    }
    return { chrToIndex: m, indexToChr: p };
  }
  async getRecordsForRange(e, t, i, r = {
    viewAsPairs: !1,
    pairAcrossChr: !1,
    maxInsertSize: 2e5
  }) {
    return j(await Te(this.streamRecordsForRange(e, t, i, r)));
  }
  async *streamRecordsForRange(e, t, i, r = {}) {
    const { signal: a } = r, o = this.chrToIndex && this.chrToIndex[e];
    let h;
    if (!(o >= 0))
      h = [];
    else if (h = await this.index.blocksForRange(o, t - 1, i, r), !h)
      throw new Error("Error in index fetch");
    for (let s = 0; s < h.length; s += 1) {
      await ee(a);
      const m = h[s].fetchedSize();
      if (m > this.chunkSizeLimit)
        throw new Error(`Too many BAM features. BAM chunk size ${m} bytes exceeds chunkSizeLimit of ${this.chunkSizeLimit}`);
    }
    const d = h.map((s) => s.fetchedSize()).reduce((s, m) => s + m, 0);
    if (d > this.fetchSizeLimit)
      throw new Error(`data size of ${d.toLocaleString()} bytes exceeded fetch size limit of ${this.fetchSizeLimit.toLocaleString()} bytes`);
    yield* this._fetchChunkFeatures(h, o, t, i, r);
  }
  async *_fetchChunkFeatures(e, t, i, r, a) {
    const { viewAsPairs: o = !1 } = a, h = [];
    let d = !1;
    for (let s = 0; s < e.length; s++) {
      const m = e[s], p = await this.featureCache.get(m.toString(), {
        chunk: m,
        opts: a
      }, a.signal), u = [];
      for (let c = 0; c < p.length; c += 1) {
        const g = p[c];
        if (g.seq_id() === t)
          if (g.get("start") >= r) {
            d = !0;
            break;
          } else
            g.get("end") >= i && u.push(g);
      }
      if (h.push(u), yield u, d)
        break;
    }
    Z(a.signal), o && (yield this.fetchPairs(t, h, a));
  }
  async fetchPairs(e, t, i) {
    const { pairAcrossChr: r = !1, maxInsertSize: a = 2e5 } = i, o = {}, h = {};
    t.map((u) => {
      const c = {};
      for (let g = 0; g < u.length; g++) {
        const _ = u[g].name(), y = u[g].id();
        c[_] || (c[_] = 0), c[_]++, h[y] = 1;
      }
      ae(c).forEach(([g, _]) => {
        _ === 1 && (o[g] = !0);
      });
    });
    const d = [];
    t.map((u) => {
      for (let c = 0; c < u.length; c++) {
        const g = u[c], _ = g.name(), y = g.get("start"), w = g._next_pos(), A = g._next_refid();
        o[_] && (r || A === e && Math.abs(y - w) < a) && d.push(this.index.blocksForRange(A, w, w + 1, i));
      }
    });
    const s = j(await Promise.all(d)).sort().filter((u, c, g) => !c || u.toString() !== g[c - 1].toString()), m = s.map((u) => u.fetchedSize()).reduce((u, c) => u + c, 0);
    if (m > this.fetchSizeLimit)
      throw new Error(`data size of ${m.toLocaleString()} bytes exceeded fetch size limit of ${this.fetchSizeLimit.toLocaleString()} bytes`);
    const p = s.map(async (u) => {
      const { data: c, cpositions: g, dpositions: _, chunk: y } = await this._readChunk({
        chunk: u,
        opts: i
      }), w = await this.readBamFeatures(c, g, _, y), A = [];
      for (let I = 0; I < w.length; I += 1) {
        const v = w[I];
        o[v.get("name")] && !h[v.id()] && A.push(v);
      }
      return A;
    });
    return j(await Promise.all(p));
  }
  async _readChunk({ chunk: e, opts: t }) {
    const i = e.fetchedSize(), { buffer: r, bytesRead: a } = await this.bam.read(O.Buffer.alloc(i), 0, i, e.minv.blockPosition, t), { buffer: o, cpositions: h, dpositions: d } = await oe(r.subarray(0, Math.min(a, i)), e);
    return { data: o, cpositions: h, dpositions: d, chunk: e };
  }
  async readBamFeatures(e, t, i, r) {
    let a = 0;
    const o = [];
    let h = 0, d = +Date.now();
    for (; a + 4 < e.length; ) {
      const s = e.readInt32LE(a), m = a + 4 + s - 1;
      if (i) {
        for (; a + r.minv.dataPosition >= i[h++]; )
          ;
        h--;
      }
      if (m < e.length) {
        const p = new Le({
          bytes: {
            byteArray: e,
            start: a,
            end: m
          },
          // the below results in an automatically calculated file-offset based ID
          // if the info for that is available, otherwise crc32 of the features
          //
          // cpositions[pos] refers to actual file offset of a bgzip block boundaries
          //
          // we multiply by (1 <<8) in order to make sure each block has a "unique"
          // address space so that data in that block could never overlap
          //
          // then the blockStart-dpositions is an uncompressed file offset from
          // that bgzip block boundary, and since the cpositions are multiplied by
          // (1 << 8) these uncompressed offsets get a unique space
          //
          // this has an extra chunk.minv.dataPosition added on because it blockStart
          // starts at 0 instead of chunk.minv.dataPosition
          //
          // the +1 is just to avoid any possible uniqueId 0 but this does not realistically happen
          fileOffset: t ? t[h] * 256 + (a - i[h]) + r.minv.dataPosition + 1 : (
            // must be slice, not subarray for buffer polyfill on web
            ce.signed(e.slice(a, m))
          )
        });
        o.push(p), this.yieldThreadTime && +Date.now() - d > this.yieldThreadTime && (await pe(1), d = +Date.now());
      }
      a = m + 1;
    }
    return o;
  }
  async hasRefSeq(e) {
    const t = this.chrToIndex && this.chrToIndex[e];
    return this.index.hasRefSeq(t);
  }
  async lineCount(e) {
    const t = this.chrToIndex && this.chrToIndex[e];
    return this.index.lineCount(t);
  }
  async indexCov(e, t, i) {
    await this.index.parse();
    const r = this.chrToIndex && this.chrToIndex[e];
    return this.index.indexCov(r, t, i);
  }
  async blocksForRange(e, t, i, r) {
    await this.index.parse();
    const a = this.chrToIndex && this.chrToIndex[e];
    return this.index.blocksForRange(a, t, i, r);
  }
}
(function() {
  (function(l) {
    (function(e) {
      var t = {
        searchParams: "URLSearchParams" in l,
        iterable: "Symbol" in l && "iterator" in Symbol,
        blob: "FileReader" in l && "Blob" in l && function() {
          try {
            return new Blob(), !0;
          } catch {
            return !1;
          }
        }(),
        formData: "FormData" in l,
        arrayBuffer: "ArrayBuffer" in l
      };
      function i(n) {
        return n && DataView.prototype.isPrototypeOf(n);
      }
      if (t.arrayBuffer)
        var r = [
          "[object Int8Array]",
          "[object Uint8Array]",
          "[object Uint8ClampedArray]",
          "[object Int16Array]",
          "[object Uint16Array]",
          "[object Int32Array]",
          "[object Uint32Array]",
          "[object Float32Array]",
          "[object Float64Array]"
        ], a = ArrayBuffer.isView || function(n) {
          return n && r.indexOf(Object.prototype.toString.call(n)) > -1;
        };
      function o(n) {
        if (typeof n != "string" && (n = String(n)), /[^a-z0-9\-#$%&'*+.^_`|~]/i.test(n))
          throw new TypeError("Invalid character in header field name");
        return n.toLowerCase();
      }
      function h(n) {
        return typeof n != "string" && (n = String(n)), n;
      }
      function d(n) {
        var f = {
          next: function() {
            var b = n.shift();
            return { done: b === void 0, value: b };
          }
        };
        return t.iterable && (f[Symbol.iterator] = function() {
          return f;
        }), f;
      }
      function s(n) {
        this.map = {}, n instanceof s ? n.forEach(function(f, b) {
          this.append(b, f);
        }, this) : Array.isArray(n) ? n.forEach(function(f) {
          this.append(f[0], f[1]);
        }, this) : n && Object.getOwnPropertyNames(n).forEach(function(f) {
          this.append(f, n[f]);
        }, this);
      }
      s.prototype.append = function(n, f) {
        n = o(n), f = h(f);
        var b = this.map[n];
        this.map[n] = b ? b + ", " + f : f;
      }, s.prototype.delete = function(n) {
        delete this.map[o(n)];
      }, s.prototype.get = function(n) {
        return n = o(n), this.has(n) ? this.map[n] : null;
      }, s.prototype.has = function(n) {
        return this.map.hasOwnProperty(o(n));
      }, s.prototype.set = function(n, f) {
        this.map[o(n)] = h(f);
      }, s.prototype.forEach = function(n, f) {
        for (var b in this.map)
          this.map.hasOwnProperty(b) && n.call(f, this.map[b], b, this);
      }, s.prototype.keys = function() {
        var n = [];
        return this.forEach(function(f, b) {
          n.push(b);
        }), d(n);
      }, s.prototype.values = function() {
        var n = [];
        return this.forEach(function(f) {
          n.push(f);
        }), d(n);
      }, s.prototype.entries = function() {
        var n = [];
        return this.forEach(function(f, b) {
          n.push([b, f]);
        }), d(n);
      }, t.iterable && (s.prototype[Symbol.iterator] = s.prototype.entries);
      function m(n) {
        if (n.bodyUsed)
          return Promise.reject(new TypeError("Already read"));
        n.bodyUsed = !0;
      }
      function p(n) {
        return new Promise(function(f, b) {
          n.onload = function() {
            f(n.result);
          }, n.onerror = function() {
            b(n.error);
          };
        });
      }
      function u(n) {
        var f = new FileReader(), b = p(f);
        return f.readAsArrayBuffer(n), b;
      }
      function c(n) {
        var f = new FileReader(), b = p(f);
        return f.readAsText(n), b;
      }
      function g(n) {
        for (var f = new Uint8Array(n), b = new Array(f.length), S = 0; S < f.length; S++)
          b[S] = String.fromCharCode(f[S]);
        return b.join("");
      }
      function _(n) {
        if (n.slice)
          return n.slice(0);
        var f = new Uint8Array(n.byteLength);
        return f.set(new Uint8Array(n)), f.buffer;
      }
      function y() {
        return this.bodyUsed = !1, this._initBody = function(n) {
          this._bodyInit = n, n ? typeof n == "string" ? this._bodyText = n : t.blob && Blob.prototype.isPrototypeOf(n) ? this._bodyBlob = n : t.formData && FormData.prototype.isPrototypeOf(n) ? this._bodyFormData = n : t.searchParams && URLSearchParams.prototype.isPrototypeOf(n) ? this._bodyText = n.toString() : t.arrayBuffer && t.blob && i(n) ? (this._bodyArrayBuffer = _(n.buffer), this._bodyInit = new Blob([this._bodyArrayBuffer])) : t.arrayBuffer && (ArrayBuffer.prototype.isPrototypeOf(n) || a(n)) ? this._bodyArrayBuffer = _(n) : this._bodyText = n = Object.prototype.toString.call(n) : this._bodyText = "", this.headers.get("content-type") || (typeof n == "string" ? this.headers.set("content-type", "text/plain;charset=UTF-8") : this._bodyBlob && this._bodyBlob.type ? this.headers.set("content-type", this._bodyBlob.type) : t.searchParams && URLSearchParams.prototype.isPrototypeOf(n) && this.headers.set("content-type", "application/x-www-form-urlencoded;charset=UTF-8"));
        }, t.blob && (this.blob = function() {
          var n = m(this);
          if (n)
            return n;
          if (this._bodyBlob)
            return Promise.resolve(this._bodyBlob);
          if (this._bodyArrayBuffer)
            return Promise.resolve(new Blob([this._bodyArrayBuffer]));
          if (this._bodyFormData)
            throw new Error("could not read FormData body as blob");
          return Promise.resolve(new Blob([this._bodyText]));
        }, this.arrayBuffer = function() {
          return this._bodyArrayBuffer ? m(this) || Promise.resolve(this._bodyArrayBuffer) : this.blob().then(u);
        }), this.text = function() {
          var n = m(this);
          if (n)
            return n;
          if (this._bodyBlob)
            return c(this._bodyBlob);
          if (this._bodyArrayBuffer)
            return Promise.resolve(g(this._bodyArrayBuffer));
          if (this._bodyFormData)
            throw new Error("could not read FormData body as text");
          return Promise.resolve(this._bodyText);
        }, t.formData && (this.formData = function() {
          return this.text().then(v);
        }), this.json = function() {
          return this.text().then(JSON.parse);
        }, this;
      }
      var w = ["DELETE", "GET", "HEAD", "OPTIONS", "POST", "PUT"];
      function A(n) {
        var f = n.toUpperCase();
        return w.indexOf(f) > -1 ? f : n;
      }
      function I(n, f) {
        f = f || {};
        var b = f.body;
        if (n instanceof I) {
          if (n.bodyUsed)
            throw new TypeError("Already read");
          this.url = n.url, this.credentials = n.credentials, f.headers || (this.headers = new s(n.headers)), this.method = n.method, this.mode = n.mode, this.signal = n.signal, !b && n._bodyInit != null && (b = n._bodyInit, n.bodyUsed = !0);
        } else
          this.url = String(n);
        if (this.credentials = f.credentials || this.credentials || "same-origin", (f.headers || !this.headers) && (this.headers = new s(f.headers)), this.method = A(f.method || this.method || "GET"), this.mode = f.mode || this.mode || null, this.signal = f.signal || this.signal, this.referrer = null, (this.method === "GET" || this.method === "HEAD") && b)
          throw new TypeError("Body not allowed for GET or HEAD requests");
        this._initBody(b);
      }
      I.prototype.clone = function() {
        return new I(this, { body: this._bodyInit });
      };
      function v(n) {
        var f = new FormData();
        return n.trim().split("&").forEach(function(b) {
          if (b) {
            var S = b.split("="), x = S.shift().replace(/\+/g, " "), E = S.join("=").replace(/\+/g, " ");
            f.append(decodeURIComponent(x), decodeURIComponent(E));
          }
        }), f;
      }
      function re(n) {
        var f = new s(), b = n.replace(/\r?\n[\t ]+/g, " ");
        return b.split(/\r?\n/).forEach(function(S) {
          var x = S.split(":"), E = x.shift().trim();
          if (E) {
            var F = x.join(":").trim();
            f.append(E, F);
          }
        }), f;
      }
      y.call(I.prototype);
      function L(n, f) {
        f || (f = {}), this.type = "default", this.status = f.status === void 0 ? 200 : f.status, this.ok = this.status >= 200 && this.status < 300, this.statusText = "statusText" in f ? f.statusText : "OK", this.headers = new s(f.headers), this.url = f.url || "", this._initBody(n);
      }
      y.call(L.prototype), L.prototype.clone = function() {
        return new L(this._bodyInit, {
          status: this.status,
          statusText: this.statusText,
          headers: new s(this.headers),
          url: this.url
        });
      }, L.error = function() {
        var n = new L(null, { status: 0, statusText: "" });
        return n.type = "error", n;
      };
      var ie = [301, 302, 303, 307, 308];
      L.redirect = function(n, f) {
        if (ie.indexOf(f) === -1)
          throw new RangeError("Invalid status code");
        return new L(null, { status: f, headers: { location: n } });
      }, e.DOMException = l.DOMException;
      try {
        new e.DOMException();
      } catch {
        e.DOMException = function(f, b) {
          this.message = f, this.name = b;
          var S = Error(f);
          this.stack = S.stack;
        }, e.DOMException.prototype = Object.create(Error.prototype), e.DOMException.prototype.constructor = e.DOMException;
      }
      function q(n, f) {
        return new Promise(function(b, S) {
          var x = new I(n, f);
          if (x.signal && x.signal.aborted)
            return S(new e.DOMException("Aborted", "AbortError"));
          var E = new XMLHttpRequest();
          function F() {
            E.abort();
          }
          E.onload = function() {
            var P = {
              status: E.status,
              statusText: E.statusText,
              headers: re(E.getAllResponseHeaders() || "")
            };
            P.url = "responseURL" in E ? E.responseURL : P.headers.get("X-Request-URL");
            var z = "response" in E ? E.response : E.responseText;
            b(new L(z, P));
          }, E.onerror = function() {
            S(new TypeError("Network request failed"));
          }, E.ontimeout = function() {
            S(new TypeError("Network request failed"));
          }, E.onabort = function() {
            S(new e.DOMException("Aborted", "AbortError"));
          }, E.open(x.method, x.url, !0), x.credentials === "include" ? E.withCredentials = !0 : x.credentials === "omit" && (E.withCredentials = !1), "responseType" in E && t.blob && (E.responseType = "blob"), x.headers.forEach(function(P, z) {
            E.setRequestHeader(z, P);
          }), x.signal && (x.signal.addEventListener("abort", F), E.onreadystatechange = function() {
            E.readyState === 4 && x.signal.removeEventListener("abort", F);
          }), E.send(typeof x._bodyInit > "u" ? null : x._bodyInit);
        });
      }
      return q.polyfill = !0, l.fetch || (l.fetch = q, l.Headers = s, l.Request = I, l.Response = L), e.Headers = s, e.Request = I, e.Response = L, e.fetch = q, Object.defineProperty(e, "__esModule", { value: !0 }), e;
    })({});
  })(typeof self < "u" ? self : he);
})();
class G {
  // uses parameter properties to automatically create fields on the class
  // https://www.typescriptlang.org/docs/handbook/classes.html#parameter-properties
  constructor(e, t, i) {
    this.record = e, this.adapter = t, this.ref = i;
  }
  _get_name() {
    return this.record.get("name");
  }
  _get_type() {
    return "match";
  }
  _get_score() {
    return this.record.get("mq");
  }
  _get_flags() {
    return this.record.flags;
  }
  _get_strand() {
    return this.record.isReverseComplemented() ? -1 : 1;
  }
  _get_pair_orientation() {
    return this.record.isPaired() ? this.record.getPairOrientation() : void 0;
  }
  _get_next_ref() {
    return this.record.isPaired() ? this.adapter.refIdToName(this.record._next_refid()) : void 0;
  }
  _get_next_pos() {
    return this.record.isPaired() ? this.record._next_pos() : void 0;
  }
  _get_next_segment_position() {
    return this.record.isPaired() ? `${this.adapter.refIdToName(this.record._next_refid())}:${this.record._next_pos() + 1}` : void 0;
  }
  _get_seq() {
    return this.record.getReadBases();
  }
  qualRaw() {
    return this.record.qualRaw();
  }
  set() {
  }
  tags() {
    const e = Object.getOwnPropertyNames(G.prototype);
    return [
      ...new Set(e.filter((t) => t.startsWith("_get_") && t !== "_get_mismatches" && t !== "_get_tags").map((t) => t.replace("_get_", "")).concat(this.record._tags()))
    ];
  }
  id() {
    return `${this.adapter.id}-${this.record.id()}`;
  }
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  get(e) {
    const t = `_get_${e}`;
    return this[t] ? this[t]() : this.record.get(e);
  }
  _get_refName() {
    return this.adapter.refIdToName(this.record.seq_id());
  }
  parent() {
  }
  children() {
  }
  pairedFeature() {
    return !1;
  }
  toJSON() {
    return {
      ...Object.fromEntries(this.tags().map((e) => [e, this.get(e)]).filter((e) => e[1] !== void 0)),
      uniqueId: this.id()
    };
  }
  _get_mismatches() {
    return fe(this.get("CIGAR"), this.get("MD"), this.get("seq"), this.ref, this.qualRaw());
  }
  _get_clipPos() {
    const e = this.get("CIGAR") || "";
    return de(e, this.get("strand"));
  }
}
class Fe extends ue.BaseFeatureDataAdapter {
  // derived classes may not use the same configuration so a custom
  // configure method allows derived classes to override this behavior
  async configurePre() {
    const e = this.getConf("bamLocation"), t = this.getConf(["index", "location"]), i = this.getConf(["index", "indexType"]), r = this.pluginManager, a = i === "CSI", o = new Be({
      bamFilehandle: N.openLocation(e, r),
      csiFilehandle: a ? N.openLocation(t, r) : void 0,
      baiFilehandle: a ? void 0 : N.openLocation(t, r),
      // chunkSizeLimit and fetchSizeLimit are more troublesome than
      // helpful, and have given overly large values on the ultra long
      // nanopore reads even with 500MB limits, so disabled with infinity
      chunkSizeLimit: 1 / 0,
      fetchSizeLimit: 1 / 0,
      yieldThreadTime: 1 / 0
    }), h = this.getConf("sequenceAdapter");
    if (h && this.getSubAdapter) {
      const { dataAdapter: d } = await this.getSubAdapter(h);
      return {
        bam: o,
        sequenceAdapter: d
      };
    } else
      return { bam: o };
  }
  async configure() {
    return this.configureP || (this.configureP = this.configurePre().catch((e) => {
      throw this.configureP = void 0, e;
    })), this.configureP;
  }
  async getHeader(e) {
    const { bam: t } = await this.configure();
    return t.getHeaderText(e);
  }
  async setupPre(e) {
    const { statusCallback: t = () => {
    } } = e || {}, { bam: i } = await this.configure();
    return this.samHeader = await V.updateStatus("Downloading index", t, async () => {
      const r = await i.getHeader(e), a = [], o = {};
      return r.filter((h) => h.tag === "SQ").forEach((h, d) => {
        h.data.forEach((s) => {
          if (s.tag === "SN") {
            const m = s.value;
            o[m] = d, a[d] = m;
          }
        });
      }), { idToName: a, nameToId: o };
    }), this.samHeader;
  }
  async setup(e) {
    return this.setupP || (this.setupP = this.setupPre(e).catch((t) => {
      throw this.setupP = void 0, t;
    })), this.setupP;
  }
  async getRefNames(e) {
    const { idToName: t } = await this.setup(e);
    return t;
  }
  async seqFetch(e, t, i) {
    const { sequenceAdapter: r } = await this.configure(), a = r;
    if (!a || !e)
      return;
    const o = a.getFeatures({
      refName: e,
      start: t,
      end: i,
      assemblyName: ""
    }), h = await le.firstValueFrom(o.pipe(me.toArray()));
    let d = "";
    if (h.sort((s, m) => s.get("start") - m.get("start")).forEach((s) => {
      const m = s.get("start"), p = s.get("end"), u = Math.max(t - m, 0), g = Math.min(i - m, p - m) - u, _ = s.get("seq") || s.get("residues");
      d += _.slice(u, u + g);
    }), d.length !== i - t)
      throw new Error(`sequence fetch failed: fetching ${e}:${(t - 1).toLocaleString()}-${i.toLocaleString()} returned ${d.length.toLocaleString()} bases, but should have returned ${(i - t).toLocaleString()}`);
    return d;
  }
  getFeatures(e, t) {
    const { refName: i, start: r, end: a, originalRefName: o } = e, { signal: h, filterBy: d, statusCallback: s = () => {
    } } = t || {};
    return ge(async (m) => {
      const { bam: p } = await this.configure();
      await this.setup(t), s("Downloading alignments");
      const u = await p.getRecordsForRange(i, r, a, t), { flagInclude: c = 0, flagExclude: g = 0, tagFilter: _, readName: y } = d || {};
      for (const w of u) {
        let A;
        w.get("MD") || (A = await this.seqFetch(o || i, w.get("start"), w.get("end")));
        const I = w.flags;
        if ((I & c) === c && !(I & g)) {
          if (_) {
            const v = w.get(_.tag);
            if (!(v === "*" ? v !== void 0 : v === _.value))
              continue;
          }
          y && w.get("name") !== y || m.next(new G(w, this, A));
        }
      }
      s(""), m.complete();
    }, h);
  }
  async getMultiRegionFeatureDensityStats(e, t) {
    const { bam: i } = await this.configure();
    if (i.index.filehandle !== "?") {
      const r = await V.bytesForRegions(e, i), a = this.getConf("fetchSizeLimit");
      return { bytes: r, fetchSizeLimit: a };
    } else
      return super.getMultiRegionFeatureDensityStats(e, t);
  }
  freeResources() {
  }
  // depends on setup being called before the BAM constructor
  refIdToName(e) {
    var t;
    return (t = this.samHeader) === null || t === void 0 ? void 0 : t.idToName[e];
  }
}
const ke = /* @__PURE__ */ Object.freeze(/* @__PURE__ */ Object.defineProperty({
  __proto__: null,
  default: Fe
}, Symbol.toStringTag, { value: "Module" }));
export {
  Be as B,
  Pe as a,
  Fe as b,
  ke as c,
  Re as p
};
