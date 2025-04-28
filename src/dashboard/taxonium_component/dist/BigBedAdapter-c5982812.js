import { B as E, P as g } from "./bbi-0b2dd278.js";
import { b8 as P, b3 as I, aX as v, aY as N, aV as _, aZ as F, l as H, aW as U, ar as V } from "./index-8ee3baa4.js";
import { Q as j } from "./index-05a6222b.js";
import { B as A, u as L } from "./util-1d4160d7.js";
import "react";
import "react-dom";
function R(p) {
  return p.filter((t) => !!t);
}
class T extends E {
  constructor() {
    super(...arguments), this.readIndicesCache = new P({
      cache: new j({ maxSize: 1 }),
      fill: async (t, e) => this._readIndices({ ...t, signal: e })
    });
  }
  readIndices(t = {}) {
    const e = "aborted" in t ? { signal: t } : t;
    return this.readIndicesCache.get(JSON.stringify(e), e, e.signal);
  }
  /*
   * retrieve unzoomed view for any scale
   * @param scale - unused
   * @param abortSignal - an optional AbortSignal to kill operation
   * @return promise for a BlockView
   */
  async getView(t, e) {
    return this.getUnzoomedView(e);
  }
  /*
   * parse the bigbed extraIndex fields
   * @param abortSignal to abort operation
   * @return a Promise for an array of Index data structure since there can be multiple extraIndexes in a bigbed, see bedToBigBed documentation
   */
  async _readIndices(t) {
    const { extHeaderOffset: e, isBigEndian: a } = await this.getHeader(t), { buffer: c } = await this.bbi.read(I.Buffer.alloc(64), 0, 64, Number(e)), o = a ? "big" : "little", l = new g().endianess(o).uint16("size").uint16("count").uint64("offset").parse(c), { count: s, offset: n } = l;
    if (s === 0)
      return [];
    const r = 20, h = r * s, { buffer: i } = await this.bbi.read(I.Buffer.alloc(h), 0, h, Number(n)), b = new g().endianess(o).int16("type").int16("fieldcount").uint64("offset").skip(4).int16("field"), f = [];
    for (let y = 0; y < s; y += 1)
      f.push(b.parse(i.subarray(y * r)));
    return f;
  }
  /*
   * perform a search in the bigbed extraIndex to find which blocks in the bigbed data to look for the
   * actual feature data
   *
   * @param name - the name to search for
   * @param opts - a SearchOptions argument with optional signal
   * @return a Promise for an array of bigbed block Loc entries
   */
  async searchExtraIndexBlocks(t, e = {}) {
    const { isBigEndian: a } = await this.getHeader(e), c = await this.readIndices(e);
    if (!c.length)
      return [];
    const o = c.map(async (l) => {
      const { offset: s, field: n } = l, { buffer: r } = await this.bbi.read(I.Buffer.alloc(32), 0, 32, Number(s), e), h = a ? "big" : "little", i = new g().endianess(h).int32("magic").int32("blockSize").int32("keySize").int32("valSize").uint64("itemCount"), { blockSize: b, keySize: f, valSize: y } = i.parse(r), z = new g().endianess(h).int8("nodeType").skip(1).int16("cnt").choice({
        tag: "nodeType",
        choices: {
          0: new g().array("leafkeys", {
            length: "cnt",
            type: new g().endianess(h).string("key", { length: f, stripNull: !0 }).uint64("offset")
          }),
          1: new g().array("keys", {
            length: "cnt",
            type: new g().endianess(h).string("key", { length: f, stripNull: !0 }).uint64("offset").uint32("length").uint32("reserved")
          })
        }
      }), k = async ($) => {
        const q = Number($), S = 4 + b * (f + y), { buffer: B } = await this.bbi.read(I.Buffer.alloc(S), 0, S, q, e), d = z.parse(B);
        if (d.leafkeys) {
          let u;
          for (let w = 0; w < d.leafkeys.length; w += 1) {
            const { key: m } = d.leafkeys[w];
            if (t.localeCompare(m) < 0 && u)
              return k(u);
            u = d.leafkeys[w].offset;
          }
          return k(u);
        }
        for (let u = 0; u < d.keys.length; u += 1)
          if (d.keys[u].key === t)
            return { ...d.keys[u], field: n };
      }, x = 32;
      return k(Number(s) + x);
    });
    return R(await Promise.all(o));
  }
  /*
   * retrieve the features from the bigbed data that were found through the lookup of the extraIndex
   * note that there can be multiple extraIndex, see the BigBed specification and the -extraIndex argument to bedToBigBed
   *
   * @param name - the name to search for
   * @param opts - a SearchOptions argument with optional signal
   * @return a Promise for an array of Feature
   */
  async searchExtraIndex(t, e = {}) {
    const a = await this.searchExtraIndexBlocks(t, e);
    if (!a.length)
      return [];
    const c = await this.getUnzoomedView(e), o = a.map((s) => new v.Observable((n) => {
      c.readFeatures(n, [s], e);
    }).pipe(N.reduce((n, r) => n.concat(r)), N.map((n) => {
      for (let r = 0; r < n.length; r += 1)
        n[r].field = s.field;
      return n;
    })));
    return (await v.firstValueFrom(v.merge(...o))).filter((s) => {
      var n;
      return ((n = s.rest) === null || n === void 0 ? void 0 : n.split("	")[(s.field || 0) - 3]) === t;
    });
  }
}
function D(p) {
  return p.get("thickStart") && p.get("blockCount") && p.get("strand") !== 0;
}
class Z extends _.BaseFeatureDataAdapter {
  async configurePre(t) {
    const e = new T({
      filehandle: F.openLocation(H.readConfObject(this.config, "bigBedLocation"), this.pluginManager)
    }), a = await e.getHeader(t), c = new A({ autoSql: a.autoSql });
    return { bigbed: e, header: a, parser: c };
  }
  async configure(t) {
    return this.cached || (this.cached = this.configurePre(t).catch((e) => {
      throw this.cached = void 0, e;
    })), this.cached;
  }
  async getRefNames(t) {
    const { header: e } = await this.configure(t);
    return Object.keys(e.refsByName);
  }
  async getHeader(t) {
    const { parser: e, header: a } = await this.configure(t), { version: c, fileType: o } = a, { fields: l, ...s } = e.autoSql;
    return {
      version: c,
      fileType: o,
      autoSql: { ...s },
      fields: Object.fromEntries(l.map(({ name: n, comment: r }) => [n, r]))
    };
  }
  getFeatures(t, e = {}) {
    const { refName: a, start: c, end: o } = t, { signal: l } = e;
    return U(async (s) => {
      try {
        const { parser: n, bigbed: r } = await this.configure(e);
        (await r.getFeatureStream(a, c, o, {
          signal: l,
          basesPerSpan: o - c
        })).pipe(N.mergeAll(), N.map((i) => {
          const b = n.parseLine(`${a}	${i.start}	${i.end}	${i.rest}`, {
            uniqueId: i.uniqueId
          }), { blockCount: f, blockSizes: y, blockStarts: z, chromStarts: k } = b;
          if (f) {
            const d = k || z || [], u = y, w = i.start;
            b.subfeatures = [];
            for (let m = 0; m < f; m += 1) {
              const C = (d[m] || 0) + w, O = C + (u[m] || 0);
              b.subfeatures.push({
                uniqueId: `${i.uniqueId}-${m}`,
                start: C,
                end: O,
                type: "block"
              });
            }
          }
          if (i.uniqueId === void 0)
            throw new Error("invalid bbi feature");
          const { chromStart: x, chromEnd: $, chrom: q, ...S } = b, B = new V({
            id: `${this.id}-${i.uniqueId}`,
            data: {
              ...S,
              start: i.start,
              end: i.end,
              refName: a
            }
          });
          return D(B) ? L(B) : B;
        })).subscribe(s);
      } catch (n) {
        s.error(n);
      }
    }, e.signal);
  }
  freeResources() {
  }
}
export {
  Z as default
};
