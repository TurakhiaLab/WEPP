import { a as b } from "./url-3d3a4b4f.js";
import { Q as w } from "./index-05a6222b.js";
import { b8 as C, aV as y, b6 as F, aW as T, t as S } from "./index-8ee3baa4.js";
import "react";
import "react-dom";
async function k(g, t, e = {}) {
  const { defaultContent: s = {} } = e;
  let i;
  try {
    return i = await t(g, { encoding: "utf8" }), JSON.parse(i);
  } catch (r) {
    if (r.code === "ENOENT" || r.status === 404)
      return s;
    throw r;
  }
}
class R {
  constructor({ readFile: t, cacheSize: e = 100 }) {
    if (this.topList = [], this.chunkCache = new C({
      cache: new w({ maxSize: e }),
      fill: this.readChunkItems.bind(this)
    }), this.readFile = t, !this.readFile)
      throw new Error('must provide a "readFile" function');
  }
  importExisting(t, e, s, i, r) {
    this.topList = t, this.attrs = e, this.start = e.makeFastGetter("Start"), this.end = e.makeFastGetter("End"), this.lazyClass = r, this.baseURL = s, this.lazyUrlTemplate = i;
  }
  binarySearch(t, e, s) {
    let i = -1, r = t.length, a;
    for (; r - i > 1; )
      a = i + r >>> 1, s(t[a]) >= e ? r = a : i = a;
    return s === this.end ? r : i;
  }
  readChunkItems(t) {
    const e = b.resolve(this.baseURL, this.lazyUrlTemplate.replace(/\{Chunk\}/gi, t));
    return k(e, this.readFile, { defaultContent: [] });
  }
  async *iterateSublist(t, e, s, i, r, a, n) {
    const o = this.attrs.makeGetter("Chunk"), u = this.attrs.makeGetter("Sublist"), l = [];
    for (let c = this.binarySearch(t, e, r); c < t.length && c >= 0 && i * a(t[c]) < i * s; c += i) {
      if (t[c][0] === this.lazyClass) {
        const d = o(t[c]), f = this.chunkCache.get(d, d).then((p) => [p, d]);
        l.push(f);
      } else
        yield [t[c], n.concat(c)];
      const h = u(t[c]);
      h && (yield* this.iterateSublist(h, e, s, i, r, a, n.concat(c)));
    }
    for (let c = 0; c < l.length; c += 1) {
      const [h, d] = await l[c];
      h && (yield* this.iterateSublist(h, e, s, i, r, a, [
        ...n,
        d
      ]));
    }
  }
  async *iterate(t, e) {
    const s = t > e ? -1 : 1, i = t > e ? this.start : this.end, r = t > e ? this.end : this.start;
    this.topList.length > 0 && (yield* this.iterateSublist(this.topList, t, e, s, i, r, [0]));
  }
  async histogram(t, e, s) {
    const i = new Array(s);
    i.fill(0);
    const r = (e - t) / s;
    for await (const a of this.iterate(t, e)) {
      const n = Math.max(0, (this.start(a) - t) / r | 0), o = Math.min(s, (this.end(a) - t) / r | 0);
      for (let u = n; u <= o; u += 1)
        i[u] += 1;
    }
    return i;
  }
}
class _ {
  constructor(t) {
    this.classes = t, this.fields = [];
    for (let e = 0; e < t.length; e += 1) {
      this.fields[e] = {};
      for (let s = 0; s < t[e].attributes.length; s += 1)
        this.fields[e][t[e].attributes[s]] = s + 1;
      t[e].proto === void 0 && (t[e].proto = {}), t[e].isArrayAttr === void 0 && (t[e].isArrayAttr = {});
    }
  }
  /**
   * @private
   */
  attrIndices(t) {
    return this.classes.map((e) => e.attributes.indexOf(t) + 1 || e.attributes.indexOf(t.toLowerCase()) + 1 || void 0);
  }
  get(t, e) {
    if (e in this.fields[t[0]])
      return t[this.fields[t[0]][e]];
    const s = e.toLowerCase();
    if (s in this.fields[t[0]])
      return t[this.fields[t[0]][s]];
    const i = this.classes[t[0]].attributes.length + 1;
    return i >= t.length || !(e in t[i]) ? e in this.classes[t[0]].proto ? this.classes[t[0]].proto[e] : void 0 : t[i][e];
  }
  makeSetter(t) {
    return (e, s) => {
      this.set(e, t, s);
    };
  }
  makeGetter(t) {
    return (e) => this.get(e, t);
  }
  makeFastGetter(t) {
    const e = this.attrIndices(t);
    return function(i) {
      if (e[i[0]] !== void 0)
        return i[e[i[0]]];
    };
  }
  // construct(self, obj, klass) {
  //   const result = new Array(self.classes[klass].length)
  //   Object.keys(obj).forEach(attr => {
  //     this.set(result, attr, obj[attr])
  //   })
  //   return result
  // }
  /**
   * Returns fast pre-compiled getter and setter functions for use with
   * Arrays that use this representation.
   * When the returned <code>get</code> and <code>set</code> functions are
   * added as methods to an Array that contains data in this
   * representation, they provide fast access by name to the data.
   *
   * @returns {Object} <code>{ get: function() {...}, set: function(val) {...} }</code>
   *
   * @example
   * var accessors = attrs.accessors();
   * var feature = get_feature_from_someplace();
   * feature.get = accessors.get;
   * // print out the feature start and end
   * console.log( feature.get('start') + ',' + feature.get('end') );
   */
  accessors() {
    return this._accessors || (this._accessors = this._makeAccessors()), this._accessors;
  }
  /**
   * @private
   */
  _makeAccessors() {
    const t = {}, e = {
      get(i) {
        const r = this.get.field_accessors[i.toLowerCase()];
        if (r)
          return r.call(this);
      },
      set(i, r) {
        const a = this.set.field_accessors[i];
        if (a)
          return a.call(this, r);
      },
      tags() {
        return s[this[0]] || [];
      }
    };
    e.get.field_accessors = {}, e.set.field_accessors = {}, this.classes.forEach((i, r) => {
      (i.attributes || []).forEach((a, n) => {
        t[a] = t[a] || [], t[a][r] = n + 1, a = a.toLowerCase(), t[a] = t[a] || [], t[a][r] = n + 1;
      });
    });
    const s = this.classes.map((i) => i.attributes);
    return Object.keys(t).forEach((i) => {
      const r = t[i];
      e.get.field_accessors[i] = r ? function() {
        return this[r[this[0]]];
      } : function() {
      };
    }), e;
  }
}
const E = _;
class L {
  constructor({ urlTemplate: t, chunkSize: e, length: s, cacheSize: i = 100, readFile: r }, a) {
    if (this.urlTemplate = t, this.chunkSize = e, this.length = s, this.baseUrl = a === void 0 ? "" : a, this.readFile = r, !r)
      throw new Error("must provide readFile callback");
    this.chunkCache = new C({
      cache: new w({ maxSize: i }),
      fill: this.getChunk.bind(this)
    });
  }
  /**
   * call the callback on one element of the array
   * @param i index
   * @param callback callback, gets called with (i, value, param)
   * @param param (optional) callback will get this as its last parameter
   */
  index(t, e, s) {
    this.range(t, t, e, void 0, s);
  }
  /**
   * async generator for the elements in the range [start,end]
   *
   * @param start index of first element to call the callback on
   * @param end index of last element to call the callback on
   */
  async *range(t, e) {
    t = Math.max(0, t), e = Math.min(e, this.length - 1);
    const s = Math.floor(t / this.chunkSize), i = Math.floor(e / this.chunkSize), r = [];
    for (let a = s; a <= i; a += 1)
      r.push(this.chunkCache.get(a, a));
    for (let a = 0; a < r.length; a += 1) {
      const [n, o] = await r[a];
      yield* this.filterChunkData(t, e, n, o);
    }
  }
  async getChunk(t) {
    let e = this.urlTemplate.replace(/\{Chunk\}/gi, t);
    this.baseUrl && (e = b.resolve(this.baseUrl, e));
    const s = await k(e, this.readFile);
    return [t, s];
  }
  *filterChunkData(t, e, s, i) {
    const r = s * this.chunkSize, a = Math.max(0, t - r), n = Math.min(e - r, this.chunkSize - 1);
    for (let o = a; o <= n; o += 1)
      yield [o + r, i[o]];
  }
}
function N() {
  return this._uniqueID;
}
function D() {
  return this._parent;
}
function x() {
  return this.get("subfeatures");
}
class z {
  constructor({ baseUrl: t, urlTemplate: e, readFile: s, cacheSize: i = 10 }) {
    if (this.baseUrl = t, this.urlTemplates = { root: e }, this.readFile = s, !this.readFile)
      throw new Error('must provide a "readFile" function argument');
    this.dataRootCache = new C({
      cache: new w({ maxSize: i }),
      fill: this.fetchDataRoot.bind(this)
    });
  }
  makeNCList() {
    return new R({ readFile: this.readFile });
  }
  loadNCList(t, e, s) {
    t.nclist.importExisting(e.intervals.nclist, t.attrs, s, e.intervals.urlTemplate, e.intervals.lazyClass);
  }
  getDataRoot(t) {
    return this.dataRootCache.get(t, t);
  }
  fetchDataRoot(t) {
    const e = b.resolve(this.baseUrl, this.urlTemplates.root.replace(/{\s*refseq\s*}/g, t));
    return k(e, this.readFile).then((s) => (
      // trackInfo = JSON.parse( trackInfo );
      this.parseTrackInfo(s, e)
    ));
  }
  parseTrackInfo(t, e) {
    const s = {
      nclist: this.makeNCList(),
      stats: {
        featureCount: t.featureCount || 0
      }
    };
    t.intervals && (s.attrs = new E(t.intervals.classes), this.loadNCList(s, t, e));
    const { histograms: i } = t;
    if (i && i.meta) {
      for (let r = 0; r < i.meta.length; r += 1)
        i.meta[r].lazyArray = new L({ ...i.meta[r].arrayParams, readFile: this.readFile }, e);
      s._histograms = i;
    }
    return s._histograms && Object.keys(s._histograms).forEach((r) => {
      s._histograms[r].forEach((n) => {
        Object.keys(n).forEach((o) => {
          typeof n[o] == "string" && String(Number(n[o])) === n[o] && (n[o] = Number(n[o]));
        });
      });
    }), s;
  }
  async getRegionStats(t) {
    return (await this.getDataRoot(t.ref)).stats;
  }
  /**
   * fetch binned counts of feature coverage in the given region.
   *
   * @param {object} query
   * @param {string} query.refName reference sequence name
   * @param {number} query.start region start
   * @param {number} query.end region end
   * @param {number} query.numBins number of bins desired in the feature counts
   * @param {number} query.basesPerBin number of bp desired in each feature counting bin
   * @returns {object} as:
   *    `{ bins: hist, stats: statEntry }`
   */
  async getRegionFeatureDensities({ refName: t, start: e, end: s, numBins: i, basesPerBin: r }) {
    const a = await this.getDataRoot(t);
    if (i)
      r = (s - e) / i;
    else if (r)
      i = Math.ceil((s - e) / r);
    else
      throw new TypeError("numBins or basesPerBin arg required for getRegionFeatureDensities");
    const o = (a._histograms.stats || []).find((h) => h.basesPerBin >= r);
    let u = a._histograms.meta[0];
    for (let h = 0; h < a._histograms.meta.length; h += 1)
      r >= a._histograms.meta[h].basesPerBin && (u = a._histograms.meta[h]);
    let l = r / u.basesPerBin;
    if (l > 0.9 && Math.abs(l - Math.round(l)) < 1e-4) {
      const h = Math.floor(e / u.basesPerBin);
      l = Math.round(l);
      const d = [];
      for (let f = 0; f < i; f += 1)
        d[f] = 0;
      for await (const [f, p] of u.lazyArray.range(h, h + l * i - 1))
        d[Math.floor((f - h) / l)] += p;
      return { bins: d, stats: o };
    }
    return { bins: await a.nclist.histogram(e, s, i), stats: o };
  }
  /**
   * Fetch features in a given region. This method is an asynchronous generator
   * yielding feature objects.
   *
   * @param {object} args
   * @param {string} args.refName reference sequence name
   * @param {number} args.start start of region. 0-based half-open.
   * @param {number} args.end end of region. 0-based half-open.
   * @yields {object}
   */
  async *getFeatures({ refName: t, start: e, end: s }) {
    const i = await this.getDataRoot(t), r = i.attrs && i.attrs.accessors();
    for await (const [a, n] of i.nclist.iterate(e, s)) {
      if (!a.decorated) {
        const o = n.join(",");
        this.decorateFeature(r, a, `${t},${o}`);
      }
      yield a;
    }
  }
  // helper method to recursively add .get and .tags methods to a feature and its
  // subfeatures
  decorateFeature(t, e, s, i) {
    e.get = t.get, e.tags = t.tags, e._uniqueID = s, e.id = N, e._parent = i, e.parent = D, e.children = x, (e.get("subfeatures") || []).forEach((r, a) => {
      this.decorateFeature(t, r, `${s}-${a}`, e);
    }), e.decorated = !0;
  }
}
const A = { refName: "seq_id" }, M = { seq_id: "refName" };
class m {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  constructor(t, e, s) {
    this.ncFeature = t, this.uniqueId = s || t.id(), this.parentHandle = e;
  }
  set() {
    throw new Error("not implemented");
  }
  jb2TagToJb1Tag(t) {
    return (A[t] || t).toLowerCase();
  }
  jb1TagToJb2Tag(t) {
    const e = t.toLowerCase();
    return M[e] || e;
  }
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  get(t) {
    const e = this.ncFeature.get(this.jb2TagToJb1Tag(t));
    return e && t === "subfeatures" ? e.map((s) => new m(s, this)) : e;
  }
  /**
   * Get an array listing which data keys are present in this feature.
   */
  tags() {
    return this.ncFeature.tags().map((t) => this.jb1TagToJb2Tag(t));
  }
  /**
   * Get the unique ID of this feature.
   */
  id() {
    return this.uniqueId;
  }
  /**
   * Get this feature's parent feature, or undefined if none.
   */
  parent() {
    return this.parentHandle;
  }
  /**
   * Get an array of child features, or undefined if none.
   */
  children() {
    return this.get("subfeatures");
  }
  toJSON() {
    const t = { uniqueId: this.id() };
    return this.ncFeature.tags().forEach((e) => {
      const s = this.jb1TagToJb2Tag(e), i = this.ncFeature.get(e);
      s === "subfeatures" ? t.subfeatures = (i || []).map((r) => new m(r, this).toJSON()) : t[s] = i;
    }), t;
  }
}
class I extends y.BaseFeatureDataAdapter {
  constructor(t, e, s) {
    super(t, e, s);
    const i = this.getConf("refNames"), r = this.getConf("rootUrlTemplate");
    this.configRefNames = i, this.nclist = new z({
      baseUrl: "",
      urlTemplate: r.uri,
      readFile: (a) => new F(String(r.baseUri ? new URL(a, r.baseUri).toString() : a)).readFile()
    });
  }
  /**
   * Fetch features for a certain region. Use getFeaturesInRegion() if you also
   * want to verify that the store has features for the given reference sequence
   * before fetching.
   * @param region -
   * @param opts - [signal] optional signalling object for aborting the fetch
   * @returns Observable of Feature objects in the region
   */
  getFeatures(t, e = {}) {
    return T(async (s) => {
      const { signal: i } = e;
      for await (const r of this.nclist.getFeatures(t, e))
        S.checkAbortSignal(i), s.next(this.wrapFeature(r));
      s.complete();
    });
  }
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  wrapFeature(t) {
    return new m(t, void 0, `${this.id}-${t.id()}`);
  }
  async hasDataForRefName(t) {
    const e = await this.nclist.getDataRoot(t);
    return !!(e && e.stats && e.stats.featureCount);
  }
  /**
   * NCList is unable to get list of ref names so returns empty
   */
  async getRefNames() {
    return this.configRefNames || [];
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
  I as default
};
