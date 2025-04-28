import { aV as o, l as f } from "./index-8ee3baa4.js";
import n from "./FromConfigAdapter-6dea6ff0.js";
import "react";
import "react-dom";
class c extends o.BaseAdapter {
  constructor(t, r, s) {
    super(t, r, s);
    const e = f.readConfObject(t, "features");
    this.features = n.makeFeatures(e || []);
  }
  /**
   * Get refName, start, and end for all features after collapsing any overlaps
   */
  async getRegions() {
    const t = [];
    for (const [r, s] of this.features) {
      let e;
      for (const a of s)
        e && e.end >= a.get("start") && e.start <= a.get("end") ? e.end = a.get("end") : (e && t.push(e), e = {
          refName: r,
          start: a.get("start"),
          end: a.get("end")
        });
      e && t.push(e);
    }
    return t.sort((r, s) => r.refName.localeCompare(s.refName)), t;
  }
  freeResources() {
  }
}
export {
  c as default
};
