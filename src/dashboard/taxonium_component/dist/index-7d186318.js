import { as as f } from "./index-8ee3baa4.js";
const S = {
  DEL: "deletion",
  INS: "insertion",
  DUP: "duplication",
  INV: "inversion",
  INVDUP: "inverted duplication",
  CNV: "copy_number_variation",
  TRA: "translocation",
  "DUP:TANDEM": "tandem_duplication",
  NON_REF: "sequence_variant",
  "*": "sequence_variant"
};
function g(n, e, s) {
  if (!e || e.length === 0)
    return ["remark", "no alternative alleles"];
  const i = /* @__PURE__ */ new Set();
  let r = /* @__PURE__ */ new Set();
  if (e.forEach((d) => {
    let [o, t] = p(n, d, s);
    o || ([o, t] = N(n, d)), o && t && (i.add(o), r.add(t));
  }), r.size > 1) {
    const d = [...r], o = new Set(d.map((t) => {
      const u = t.split("->");
      return u[1] ? u[0] : t;
    }));
    r = new Set([...o].map((t) => {
      const u = d.map((c) => {
        const l = c.split("-> ");
        return l[1] && l[0] === t ? l[1] : "";
      }).filter((c) => !!c);
      return u.length ? t + "-> " + u.join(",") : t;
    }));
  }
  return i.size ? [[...i].join(","), [...r].join(",")] : [];
}
function p(n, e, s) {
  if (typeof e == "string" && !e.startsWith("<"))
    return [];
  let i = S[e];
  if (!i && s.getMetadata("ALT", e) && (i = "sequence_variant"), i)
    return [i, e];
  const r = e.split(":");
  return r.length > 1 ? p(n, `<${r.slice(0, -1).join(":")}>`, s) : [];
}
function N(n, e) {
  return f(e) ? ["breakend", e] : n.length === 1 && e.length === 1 ? ["SNV", a("SNV", n, e)] : e === "<INS>" ? ["insertion", e] : e === "<DEL>" ? ["deletion", e] : e === "<INV>" ? ["deletion", e] : e === "<TRA>" ? ["translocation", e] : e.includes("<") ? ["sv", e] : n.length === e.length ? n.split("").reverse().join("") === e ? ["inversion", a("inversion", n, e)] : ["substitution", a("substitution", n, e)] : n.length <= e.length ? ["insertion", a("insertion", n, e)] : n.length > e.length ? ["deletion", a("deletion", n, e)] : ["indel", a("indel", n, e)];
}
function a(n, e, s) {
  return `${n} ${e} -> ${s}`;
}
class A {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  constructor(e) {
    this.variant = e.variant, this.parser = e.parser, this.data = this.dataFromVariant(this.variant), this._id = e.id;
  }
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  get(e) {
    return e === "samples" ? this.variant.SAMPLES : this.data[e] || this.variant[e];
  }
  set() {
  }
  parent() {
  }
  children() {
  }
  tags() {
    return [...Object.keys(this.data), ...Object.keys(this.variant), "samples"];
  }
  id() {
    return this._id;
  }
  dataFromVariant(e) {
    const { REF: s, ALT: i, POS: r, CHROM: d, INFO: o, ID: t } = e, u = r - 1, [c, l] = g(s, i, this.parser), m = i == null ? void 0 : i.some((h) => h === "<TRA>"), v = i == null ? void 0 : i.some((h) => h.includes("<"));
    return {
      refName: d,
      start: u,
      end: v && o.END && !m ? +o.END[0] : u + s.length,
      description: l,
      type: c,
      name: t == null ? void 0 : t.join(","),
      aliases: t && t.length > 1 ? e.ID.slice(1) : void 0
    };
  }
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  toJSON() {
    return {
      uniqueId: this._id,
      ...this.variant,
      ...this.data,
      samples: this.variant.SAMPLES
    };
  }
}
export {
  A as V
};
