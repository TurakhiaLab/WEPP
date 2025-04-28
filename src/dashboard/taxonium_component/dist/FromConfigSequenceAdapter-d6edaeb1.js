import { aW as n, aX as o, aY as f, ar as u } from "./index-8ee3baa4.js";
import d from "./FromConfigAdapter-6dea6ff0.js";
import "react";
import "react-dom";
class l extends d {
  /**
   * Fetch features for a certain region
   * @param region - Region
   * @returns Observable of Feature objects in the region
   */
  getFeatures(e) {
    return n(async (s) => {
      const t = (await o.firstValueFrom(super.getFeatures(e).pipe(f.toArray())))[0];
      s.next(new u({
        ...t.toJSON(),
        uniqueId: `${t.id()}:${e.start}-${e.end}`,
        end: e.end,
        start: e.start,
        seq: t.get("seq").slice(Math.max(e.start - t.get("start"), 0), Math.max(e.end - t.get("start"), 0))
      })), s.complete();
    });
  }
  /**
   * Get refName, start, and end for all features after collapsing any overlaps
   */
  async getRegions() {
    const e = [];
    for (const [s, r] of this.features) {
      let t;
      for (const a of r)
        t && t.end >= a.get("start") && t.start <= a.get("end") ? t.end = a.get("end") : (t && e.push(t), t = {
          refName: s,
          start: a.get("start"),
          end: a.get("end")
        });
      t && e.push(t);
    }
    return e;
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
  l as default
};
