import n from "react";
import { o as C, q as f, u as p, b, T as s, l as m, s as h, C as y, c as E, B as k } from "./index-8ee3baa4.js";
import "react-dom";
function v(e) {
  return e.length > 20 ? `${e.slice(0, 20)}...` : e;
}
const D = f()((e) => ({
  connectionContainer: {
    width: 500,
    margin: e.spacing(4)
  }
}));
function T({ session: e, handleClose: a, breakConnection: u }) {
  const { classes: g } = D(), { connections: c, connectionInstances: d = [] } = e;
  return n.createElement(
    p.Dialog,
    { open: !0, onClose: a, maxWidth: "lg", title: "Turn on/off connections" },
    n.createElement(
      b,
      null,
      n.createElement(s, null, "Use the checkbox to turn on/off connections"),
      n.createElement(
        "div",
        { className: g.connectionContainer },
        c.map((o) => {
          const l = m.readConfObject(o, "name"), i = m.readConfObject(o, "assemblyNames"), r = d.find((t) => l === t.name);
          return n.createElement(h, { key: o.connectionId, control: n.createElement(y, { checked: !!r, onChange: () => {
            var t;
            r ? u(o) : (t = e.makeConnection) === null || t === void 0 || t.call(e, o);
          }, color: "primary" }), label: `${l} ${i.length ? "(" + v(i.join(",")) + ")" : ""}` });
        }),
        c.length ? null : n.createElement(s, null, "No connections found")
      )
    ),
    n.createElement(
      E,
      null,
      n.createElement(k, { onClick: () => a(), variant: "contained", color: "primary" }, "Close")
    )
  );
}
const _ = C(T);
export {
  _ as default,
  v as ellipses
};
