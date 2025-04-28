import { aV as d, t as l, aW as u, aX as p, aY as g } from "./index-8ee3baa4.js";
import "react";
import "react-dom";
function m(n) {
  const a = n.slice(n.lastIndexOf("/") + 1);
  return a.slice(0, a.lastIndexOf("."));
}
class A extends d.BaseFeatureDataAdapter {
  async getAdapters() {
    const a = this.getSubAdapter;
    if (!a)
      throw new Error("no getSubAdapter available");
    let e = this.getConf("subadapters");
    return e != null && e.length || (e = this.getConf("bigWigs").map((s) => ({
      type: "BigWigAdapter",
      source: m(s),
      bigWigLocation: {
        uri: s
      }
    }))), Promise.all(
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      e.map(async (t) => {
        const s = (await a(t)).dataAdapter;
        return {
          source: t.name || s.id,
          ...t,
          dataAdapter: s
        };
      })
    );
  }
  // note: can't really have dis-agreeing refNames
  async getRefNames(a) {
    const e = await this.getAdapters(), t = await Promise.all(e.map((s) => s.dataAdapter.getRefNames(a)));
    return [...new Set(t.flat())];
  }
  async getGlobalStats(a) {
    const e = await this.getAdapters(), t = (await Promise.all(
      // @ts-expect-error
      e.map((r) => {
        var c, o;
        return (o = (c = r.dataAdapter).getGlobalStats) === null || o === void 0 ? void 0 : o.call(c, a);
      })
    )).filter((r) => !!r), s = l.min(t.map((r) => r.scoreMin)), i = l.max(t.map((r) => r.scoreMax));
    return { scoreMin: s, scoreMax: i };
  }
  getFeatures(a, e = {}) {
    return u(async (t) => {
      const s = await this.getAdapters();
      p.merge(...s.map((i) => i.dataAdapter.getFeatures(a, e).pipe(g.map((r) => (
        // add source field if it does not exist
        r.get("source") ? r : new l.SimpleFeature({
          ...r.toJSON(),
          uniqueId: `${i.source}-${r.id()}`,
          source: i.source
        })
      ))))).subscribe(t);
    }, e.signal);
  }
  // always render bigwig instead of calculating a feature density for it
  async getMultiRegionFeatureDensityStats(a) {
    return { featureDensity: 0 };
  }
  // in another adapter type, this could be dynamic depending on region or
  // something, but it is static for this particular multi-wiggle adapter type
  async getSources() {
    return (await this.getAdapters()).map(({ dataAdapter: e, source: t, name: s, ...i }) => ({
      name: t,
      __name: s,
      ...i
    }));
  }
  freeResources() {
  }
}
A.capabilities = [
  "hasResolution",
  "hasLocalStats",
  "hasGlobalStats"
];
export {
  A as default
};
