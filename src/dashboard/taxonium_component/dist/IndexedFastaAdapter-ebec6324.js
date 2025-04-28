import { b5 as d, b3 as x, aV as b, b8 as _, bg as I, aZ as l, aW as L, t as q } from "./index-8ee3baa4.js";
function y(r, e) {
  return r.offset + r.lineBytes * Math.floor(e / r.lineLength) + e % r.lineLength;
}
async function F(r, e) {
  const t = await r.readFile(e);
  if (!(t && t.length))
    throw new Error("No data read from FASTA index (FAI) file");
  let a = 0, n;
  const s = t.toString("utf8").split(/\r?\n/).filter((i) => /\S/.test(i)).map((i) => i.split("	")).filter((i) => i[0] !== "").map((i) => ((!n || n.name !== i[0]) && (n = { name: i[0], id: a }, a += 1), {
    id: n.id,
    name: i[0],
    length: +i[1],
    start: 0,
    end: +i[1],
    offset: +i[2],
    lineLength: +i[3],
    lineBytes: +i[4]
  }));
  return {
    name: Object.fromEntries(s.map((i) => [i.name, i])),
    id: Object.fromEntries(s.map((i) => [i.id, i]))
  };
}
let z = class {
  constructor({ fasta: e, fai: t, path: a, faiPath: n, chunkSizeLimit: s = 1e6 }) {
    if (e)
      this.fasta = e;
    else if (a)
      this.fasta = new d(a);
    else
      throw new Error("Need to pass filehandle for fasta or path to localfile");
    if (t)
      this.fai = t;
    else if (n)
      this.fai = new d(n);
    else if (a)
      this.fai = new d(`${a}.fai`);
    else
      throw new Error("Need to pass filehandle for  or path to localfile");
    this.chunkSizeLimit = s;
  }
  async _getIndexes(e) {
    return this.indexes || (this.indexes = F(this.fai, e)), this.indexes;
  }
  /**
   * @returns {array[string]} array of string sequence
   * names that are present in the index, in which the
   * array index indicates the sequence ID, and the value
   * is the sequence name
   */
  async getSequenceNames(e) {
    return Object.keys((await this._getIndexes(e)).name);
  }
  /**
   * @returns {array[string]} array of string sequence
   * names that are present in the index, in which the
   * array index indicates the sequence ID, and the value
   * is the sequence name
   */
  async getSequenceSizes(e) {
    const t = {}, a = await this._getIndexes(e), n = Object.values(a.id);
    for (let s = 0; s < n.length; s += 1)
      t[n[s].name] = n[s].length;
    return t;
  }
  /**
   * @returns {array[string]} array of string sequence
   * names that are present in the index, in which the
   * array index indicates the sequence ID, and the value
   * is the sequence name
   */
  async getSequenceSize(e, t) {
    var a;
    return (a = (await this._getIndexes(t)).name[e]) === null || a === void 0 ? void 0 : a.length;
  }
  /**
   *
   * @param {string} name
   * @returns {Promise[boolean]} true if the file contains the given reference sequence name
   */
  async hasReferenceSequence(e, t) {
    return !!(await this._getIndexes(t)).name[e];
  }
  /**
   *
   * @param {number} seqId
   * @param {number} min
   * @param {number} max
   */
  async getResiduesById(e, t, a, n) {
    const s = (await this._getIndexes(n)).id[e];
    if (s)
      return this._fetchFromIndexEntry(s, t, a, n);
  }
  /**
   * @param {string} seqName
   * @param {number} min
   * @param {number} max
   */
  async getResiduesByName(e, t, a, n) {
    const s = (await this._getIndexes(n)).name[e];
    if (s)
      return this._fetchFromIndexEntry(s, t, a, n);
  }
  //alias for getResiduesByName
  async getSequence(e, t, a, n) {
    return this.getResiduesByName(e, t, a, n);
  }
  async _fetchFromIndexEntry(e, t = 0, a, n) {
    let s = a;
    if (t < 0)
      throw new TypeError("regionStart cannot be less than 0");
    if ((s === void 0 || s > e.length) && (s = e.length), t >= s)
      return "";
    const i = y(e, t), o = y(e, s) - i;
    if (o > this.chunkSizeLimit)
      throw new Error(`data size of ${o.toLocaleString()} bytes exceeded chunk size limit of ${this.chunkSizeLimit.toLocaleString()} bytes`);
    const c = x.Buffer.allocUnsafe(o);
    return await this.fasta.read(c, 0, o, i, n), c.toString("utf8").replace(/\s+/g, "");
  }
};
class E extends b.BaseSequenceAdapter {
  constructor() {
    super(...arguments), this.seqCache = new _({
      cache: new I({ maxSize: 200 }),
      fill: async (e, t) => {
        const { refName: a, start: n, end: s, fasta: i } = e;
        return i.getSequence(a, n, s, { ...e, signal: t });
      }
    });
  }
  async getRefNames(e) {
    const { fasta: t } = await this.setup();
    return t.getSequenceNames(e);
  }
  async getRegions(e) {
    const { fasta: t } = await this.setup(), a = await t.getSequenceSizes(e);
    return Object.keys(a).map((n) => ({
      refName: n,
      start: 0,
      end: a[n]
    }));
  }
  async setupPre() {
    const e = this.getConf("fastaLocation"), t = this.getConf("faiLocation");
    return {
      fasta: new z({
        fasta: l.openLocation(e, this.pluginManager),
        fai: l.openLocation(t, this.pluginManager)
      })
    };
  }
  async getHeader() {
    const e = this.getConf("metadataLocation");
    return e.uri === "" || e.uri === "/path/to/fa.metadata.yaml" ? null : l.openLocation(e, this.pluginManager).readFile("utf8");
  }
  async setup() {
    return this.setupP || (this.setupP = this.setupPre().catch((e) => {
      throw this.setupP = void 0, e;
    })), this.setupP;
  }
  getFeatures(e, t) {
    const { refName: a, start: n, end: s } = e;
    return L(async (i) => {
      const { fasta: o } = await this.setup(), c = await o.getSequenceSize(a, t), h = c !== void 0 ? Math.min(c, s) : s, g = [], u = 128e3, m = n - n % u, w = s + (u - s % u);
      for (let f = m; f < w; f += u) {
        const S = {
          refName: a,
          start: f,
          end: f + u
        };
        g.push(this.seqCache.get(JSON.stringify(S), { ...S, fasta: o }, t == null ? void 0 : t.signal));
      }
      const p = (await Promise.all(g)).join("").slice(n - m).slice(0, s - n);
      p && i.next(new q.SimpleFeature({
        id: `${a} ${n}-${h}`,
        data: { refName: a, start: n, end: h, seq: p }
      })), i.complete();
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
const j = /* @__PURE__ */ Object.freeze(/* @__PURE__ */ Object.defineProperty({
  __proto__: null,
  default: E
}, Symbol.toStringTag, { value: "Module" }));
export {
  z as I,
  E as a,
  j as b
};
