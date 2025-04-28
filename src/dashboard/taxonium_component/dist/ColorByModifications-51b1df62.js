import e from "react";
import { o as m, q as d, u, b as y, T as a, x as f, c as p, B as l } from "./index-8ee3baa4.js";
import "react-dom";
const g = d()((n) => ({
  table: {
    border: "1px solid #888",
    margin: n.spacing(4),
    "& td": {
      padding: n.spacing(1)
    }
  }
}));
function r({ modifications: n }) {
  const { classes: i } = g();
  return e.createElement(
    "table",
    { className: i.table },
    e.createElement("tbody", null, n.map(([o, t]) => e.createElement(
      "tr",
      { key: o },
      e.createElement("td", null, o),
      e.createElement("td", null, t),
      e.createElement("td", { style: {
        width: "1em",
        background: t
      } })
    )))
  );
}
function h(n) {
  const { model: i, handleClose: o } = n, { colorBy: t, modificationTagMap: c } = i, s = [...c.entries()];
  return e.createElement(
    u.Dialog,
    { open: !0, onClose: o, title: "Color by modifications" },
    e.createElement(
      y,
      null,
      e.createElement(a, null, "You can choose to color the modifications in the BAM/CRAM MM/ML specification using this dialog. Choosing modifications colors the modified positions and can color multiple modification types. Choosing the methylation setting colors methylated and unmethylated CpG."),
      e.createElement(a, null, "Note: you can revisit this dialog to see the current mapping of colors to modification type for the modification coloring mode"),
      e.createElement(
        "div",
        { style: { margin: 20 } },
        (t == null ? void 0 : t.type) === "modifications" ? e.createElement("div", null, s.length ? e.createElement(
          e.Fragment,
          null,
          "Current modification-type-to-color mapping",
          e.createElement(r, { modifications: [...c.entries()] })
        ) : e.createElement(
          "div",
          null,
          e.createElement(a, null, "Note: color by modifications is already enabled. Loading current modifications..."),
          e.createElement(f, { size: 15 })
        )) : null,
        (t == null ? void 0 : t.type) === "methylation" ? e.createElement(r, { modifications: [
          ["methylated", "red"],
          ["unmethylated", "blue"]
        ] }) : null
      ),
      e.createElement(
        p,
        null,
        e.createElement(l, { variant: "contained", onClick: () => {
          i.setColorScheme({ type: "modifications" }), o();
        } }, "Modifications"),
        e.createElement(l, { variant: "contained", onClick: () => {
          i.setColorScheme({ type: "methylation" }), o();
        } }, "Methylation"),
        e.createElement(l, { variant: "contained", color: "secondary", onClick: () => o() }, "Cancel")
      )
    )
  );
}
const v = m(h);
export {
  v as default
};
