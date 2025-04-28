import e, { useState as d } from "react";
import { o as B, q as P, u as R, b as S, T as u, Z as A, P as p, n as g, c as _, B as F } from "./index-8ee3baa4.js";
import "react-dom";
const D = P()((n) => ({
  paper: {
    padding: n.spacing(2),
    margin: n.spacing(2)
  },
  field: {
    margin: n.spacing(2)
  }
})), L = [
  "read paired",
  "read mapped in proper pair",
  "read unmapped",
  "mate unmapped",
  "read reverse strand",
  "mate reverse strand",
  "first in pair",
  "second in pair",
  "not primary alignment",
  "read fails platform/vendor quality checks",
  "read is PCR or optical duplicate",
  "supplementary alignment"
];
function b(n) {
  const { flag: l = 0, setFlag: r } = n;
  return e.createElement(
    e.Fragment,
    null,
    e.createElement(g, { type: "number", value: l, onChange: (o) => r(+o.target.value) }),
    L.map((o, i) => {
      const t = l & 1 << i, a = `${o}_${t}`;
      return e.createElement(
        "div",
        { key: a },
        e.createElement("input", { type: "checkbox", checked: !!t, onChange: (m) => {
          m.target.checked ? r(l | 1 << i) : r(l & ~(1 << i));
        } }),
        e.createElement("label", { htmlFor: a }, o)
      );
    })
  );
}
function I(n) {
  var l, r;
  const { model: o, handleClose: i } = n, { classes: t } = D(), { filterBy: a } = o, [m, N] = d(a == null ? void 0 : a.flagInclude), [v, k] = d(a == null ? void 0 : a.flagExclude), [s, C] = d(((l = a == null ? void 0 : a.tagFilter) === null || l === void 0 ? void 0 : l.tag) || ""), [f, x] = d(((r = a == null ? void 0 : a.tagFilter) === null || r === void 0 ? void 0 : r.value) || ""), [E, T] = d((a == null ? void 0 : a.readName) || ""), h = s.match(/^[A-Za-z][A-Za-z0-9]$/), y = "https://broadinstitute.github.io/picard/explain-flags.html";
  return e.createElement(
    R.Dialog,
    { open: !0, onClose: i, title: "Filter options" },
    e.createElement(
      S,
      null,
      e.createElement(
        u,
        null,
        "Set filter bitmask options. Refer to ",
        e.createElement(A, { href: y }, y),
        " ",
        "for details"
      ),
      e.createElement(
        p,
        { className: t.paper, variant: "outlined" },
        e.createElement(
          "div",
          { style: { display: "flex" } },
          e.createElement(
            "div",
            null,
            e.createElement(u, null, "Read must have ALL these flags"),
            e.createElement(b, { flag: m, setFlag: N })
          ),
          e.createElement(
            "div",
            null,
            e.createElement(u, null, "Read must have NONE of these flags"),
            e.createElement(b, { flag: v, setFlag: k })
          )
        )
      ),
      e.createElement(
        p,
        { className: t.paper, variant: "outlined" },
        e.createElement(u, null, "Filter by tag name and value. Use * in the value field to get all reads containing any value for that tag. Example: filter tag name SA with value * to get all split/supplementary reads. Other examples include HP for haplotype, or RG for read group"),
        e.createElement(g, { className: t.field, value: s, onChange: (c) => C(c.target.value), placeholder: "Enter tag name", inputProps: {
          maxLength: 2,
          "data-testid": "color-tag-name-input"
        }, error: s.length === 2 && !h, helperText: s.length === 2 && !h ? "Not a valid tag" : "", "data-testid": "color-tag-name" }),
        e.createElement(g, { className: t.field, value: f, onChange: (c) => x(c.target.value), placeholder: "Enter tag value", inputProps: {
          "data-testid": "color-tag-name-input"
        }, "data-testid": "color-tag-value" })
      ),
      e.createElement(
        p,
        { className: t.paper, variant: "outlined" },
        e.createElement(u, null, "Filter by read name"),
        e.createElement(g, { className: t.field, value: E, onChange: (c) => T(c.target.value), placeholder: "Enter read name", inputProps: {
          "data-testid": "color-tag-readname-input"
        }, "data-testid": "color-tag-readname" })
      ),
      e.createElement(
        _,
        null,
        e.createElement(F, { variant: "contained", color: "primary", autoFocus: !0, type: "submit", onClick: () => {
          o.setFilterBy({
            flagInclude: m,
            flagExclude: v,
            readName: E,
            tagFilter: s !== "" ? {
              tag: s,
              value: f
            } : void 0
          }), i();
        } }, "Submit"),
        e.createElement(F, { variant: "contained", color: "secondary", onClick: () => i() }, "Cancel")
      )
    )
  );
}
const z = B(I);
export {
  z as default
};
