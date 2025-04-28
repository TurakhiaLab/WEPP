import "./long-a46ae85c.js";
import { b3 as f, b9 as g } from "./index-8ee3baa4.js";
import { B as k, a as I, p as A, b as B } from "./BamAdapter-ebc7a99d.js";
import "react";
import "react-dom";
async function w(h, t) {
  const s = await Promise.all(h.map(async (e) => {
    const { url: a, headers: n } = e;
    if (a.startsWith("data:"))
      return f.Buffer.from(a.split(",")[1], "base64");
    {
      const { referer: c, ...o } = n, r = await fetch(a, {
        ...t,
        headers: { ...t.headers, ...o }
      });
      if (!r.ok)
        throw new Error(`Failed to fetch ${r.statusText}`);
      return f.Buffer.from(await r.arrayBuffer());
    }
  }));
  return f.Buffer.concat(await Promise.all(s.map((e) => g(e))));
}
class T extends k {
  constructor(t) {
    super({ bamFilehandle: "?", baiFilehandle: "?" }), this.baseUrl = t.baseUrl, this.trackId = t.trackId;
  }
  async *streamRecordsForRange(t, s, e, a = {
    viewAsPairs: !1,
    pairAcrossChr: !1,
    maxInsertSize: 2e5
  }) {
    const c = `${`${this.baseUrl}/${this.trackId}`}?referenceName=${t}&start=${s}&end=${e}&format=BAM`, o = this.chrToIndex && this.chrToIndex[t], r = await fetch(c, { ...a });
    if (!r.ok)
      throw new Error(r.statusText);
    const i = await r.json(), l = {
      buffer: await w(i.htsget.urls.slice(1), a),
      chunk: { minv: { dataPosition: 0 } },
      toString() {
        return `${t}_${s}_${e}`;
      }
    };
    yield* this._fetchChunkFeatures(
      // @ts-ignore
      [l],
      o,
      s,
      e,
      a
    );
  }
  //@ts-ignore
  async _readChunk(t) {
    const { chunk: s } = t, { buffer: e, chunk: a } = s;
    return { data: e, cpositions: null, dpositions: null, chunk: a };
  }
  async getHeader(t = {}) {
    const s = `${this.baseUrl}/${this.trackId}?referenceName=na&class=header`, e = await fetch(s, t);
    if (!e.ok)
      throw new Error(`Failed to fetch ${e.statusText}`);
    const a = await e.json(), n = await w(a.htsget.urls, t);
    if (n.readInt32LE(0) !== I)
      throw new Error("Not a BAM file");
    const c = n.readInt32LE(4), o = n.toString("utf8", 8, 8 + c), r = A(o), i = [], d = {};
    return r.filter((u) => u.tag === "SQ").forEach((u, m) => {
      u.data.forEach((b) => {
        if (b.tag === "SN") {
          const p = b.value;
          d[p] = m, i[m] = p;
        }
      });
    }), this.chrToIndex = d, this.indexToChr = i, r;
  }
}
class F extends B {
  async configurePre() {
    const t = this.getConf("htsgetBase"), s = this.getConf("htsgetTrackId"), e = new T({
      baseUrl: t,
      trackId: s
    }), a = this.getConf("sequenceAdapter");
    if (a && this.getSubAdapter) {
      const n = await this.getSubAdapter(a);
      return {
        bam: e,
        sequenceAdapter: n.dataAdapter
      };
    }
    return { bam: e };
  }
}
export {
  F as default
};
