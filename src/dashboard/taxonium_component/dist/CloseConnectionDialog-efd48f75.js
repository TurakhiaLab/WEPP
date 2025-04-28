import e from "react";
import { o as m, D as u, a as C, b as E, aE as c, a3 as g, a4 as p, c as D, B as r } from "./index-8ee3baa4.js";
import "react-dom";
function y({ modalInfo: t = {}, setModalInfo: l }) {
  const { name: i, dereferenceTypeCount: o, safelyBreakConnection: n } = t;
  return e.createElement(
    u,
    { open: !0 },
    e.createElement(
      C,
      null,
      'Close connection "',
      i,
      '"'
    ),
    e.createElement(
      E,
      null,
      o ? e.createElement(
        e.Fragment,
        null,
        e.createElement(c, null, "Closing this connection will close:"),
        e.createElement(g, null, Object.entries(o).map(([a, s]) => e.createElement(p, { key: a }, `${s} ${a}`)))
      ) : null,
      e.createElement(c, null, "Are you sure you want to close this connection?")
    ),
    e.createElement(
      D,
      null,
      e.createElement(r, { onClick: () => l(), color: "primary" }, "Cancel"),
      e.createElement(r, { variant: "contained", onClick: t ? () => {
        n == null || n(), l();
      } : () => {
      }, color: "primary" }, "OK")
    )
  );
}
const b = m(y);
export {
  b as default
};
