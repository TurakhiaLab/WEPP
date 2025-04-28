import { B as b, a as w } from "./bbi-0b2dd278.js";
import { aV as y, aZ as B, t as P, bi as S, aW as x } from "./index-8ee3baa4.js";
import "./index-05a6222b.js";
import "react";
import "react-dom";
class N extends b {
  /**
   * Retrieves a BlockView of a specific zoomLevel
   *
   * @param scale - number
   * @param opts - An object containing basesPerSpan (e.g. pixels per basepair) or scale used to infer the zoomLevel to use
   */
  async getView(t, e) {
    const { zoomLevels: a, refsByName: r, fileSize: u, isBigEndian: c, uncompressBufSize: l } = await this.getHeader(e), d = 1 / t;
    let o = a.length;
    u || (o -= 1);
    for (let n = o; n >= 0; n -= 1) {
      const s = a[n];
      if (s && s.reductionLevel <= 2 * d) {
        const g = Number(s.indexOffset);
        return new w(this.bbi, r, g, c, l > 0, "summary");
      }
    }
    return this.getUnzoomedView(e);
  }
}
class v extends y.BaseFeatureDataAdapter {
  async setupPre(t) {
    const { statusCallback: e = () => {
    } } = t || {}, a = new N({
      filehandle: B.openLocation(this.getConf("bigWigLocation"), this.pluginManager)
    }), r = await P.updateStatus("Downloading bigwig header", e, () => a.getHeader(t));
    return { bigwig: a, header: r };
  }
  async setup(t) {
    return this.setupP || (this.setupP = this.setupPre(t).catch((e) => {
      throw this.setupP = void 0, e;
    })), this.setupP;
  }
  async getRefNames(t) {
    const { header: e } = await this.setup(t);
    return Object.keys(e.refsByName);
  }
  async refIdToName(t) {
    var e;
    const { header: a } = await this.setup();
    return (e = a.refsByNumber[t]) === null || e === void 0 ? void 0 : e.name;
  }
  async getGlobalStats(t) {
    const { header: e } = await this.setup(t);
    return S(e.totalSummary);
  }
  getFeatures(t, e = {}) {
    const { refName: a, start: r, end: u } = t, { bpPerPx: c = 0, signal: l, resolution: d = 1, statusCallback: o = () => {
    } } = e;
    return x(async (n) => {
      o("Downloading bigwig data");
      const s = this.getConf("source"), { bigwig: g } = await this.setup(e), m = await g.getFeatures(a, r, u, {
        ...e,
        basesPerSpan: c / d
      });
      for (const i of m) {
        s && (i.source = s);
        const f = `${s}:${t.refName}:${i.start}-${i.end}`;
        i.refName = a, i.uniqueId = f, n.next({
          // @ts-expect-error
          get: (p) => i[p],
          id: () => f,
          // @ts-expect-error
          toJSON: () => i
        });
      }
      n.complete();
    }, l);
  }
  // always render bigwig instead of calculating a feature density for it
  async getMultiRegionFeatureDensityStats(t) {
    return { featureDensity: 0 };
  }
  freeResources() {
  }
}
v.capabilities = [
  "hasResolution",
  "hasLocalStats",
  "hasGlobalStats"
];
export {
  v as default
};
