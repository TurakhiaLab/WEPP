import e from "react";
import { o as g, q as E, u as p, b as f, T as a, l as C, I as r, U as s, a9 as y, c as D, B as b } from "./index-8ee3baa4.js";
import "react-dom";
const h = E()((n) => ({
  connectionContainer: {
    margin: n.spacing(4),
    width: 500
  }
}));
function k({ session: n, handleClose: l, breakConnection: m }) {
  const { classes: u } = h(), { adminMode: d, connections: c, sessionConnections: t } = n;
  return e.createElement(
    p.Dialog,
    { open: !0, onClose: l, maxWidth: "lg", title: "Delete connections" },
    e.createElement(
      f,
      null,
      e.createElement(a, null, "Click the X icon to delete the connection from your config completely"),
      e.createElement(
        "div",
        { className: u.connectionContainer },
        c.map((o) => {
          const i = C.readConfObject(o, "name");
          return e.createElement(
            a,
            { key: `conn-${i}` },
            d || t != null && t.includes(o) ? e.createElement(
              r,
              { onClick: () => m(o, !0) },
              e.createElement(s, { color: "error" })
            ) : e.createElement(
              y,
              { title: "Unable to delete connection in config file as non-admin user" },
              e.createElement(
                r,
                null,
                e.createElement(s, { color: "disabled" })
              )
            ),
            i
          );
        }),
        c.length ? null : e.createElement(a, null, "No connections found")
      )
    ),
    e.createElement(
      D,
      null,
      e.createElement(b, { onClick: () => l(), variant: "contained", color: "primary" }, "Close")
    )
  );
}
const T = g(k);
export {
  T as default
};
