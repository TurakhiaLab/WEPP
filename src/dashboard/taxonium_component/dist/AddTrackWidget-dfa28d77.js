import e, { useState as m, useEffect as L } from "react";
import { q as v, o as k, P as _, a2 as B, ab as H, ad as V, a3 as z, a4 as W, n as b, a5 as F, I as $, m as K, p as Y, t as d, M as D, at as Z, T as h, Z as w, au as Q, u as N, ae as U, s as X, C as ee, H as j, av as te, aw as ae, ax as ne, ay as re, B as I, az as se, l as ce, aA as le, a6 as oe } from "./index-8ee3baa4.js";
import "react-dom";
const ie = v()((t) => ({
  paper: {
    display: "flex",
    flexDirection: "column",
    padding: t.spacing(1)
  },
  card: {
    marginTop: t.spacing(1)
  }
})), ue = k(function({ model: t }) {
  const { classes: a } = ie(), [r, c] = m(""), [s, l] = m(""), [o, n] = m(["Name", "ID"]), [i, u] = m(["CDS", "exon"]), E = [
    {
      label: "Indexing attributes",
      values: o
    },
    {
      label: "Feature types to exclude",
      values: i
    }
  ];
  return L(() => {
    t.setTextIndexingConf({ attributes: o, exclude: i });
  }, [t, o, i]), e.createElement(
    _,
    { className: a.paper },
    e.createElement(B, null, "Indexing configuration"),
    E.map((g, p) => e.createElement(
      H,
      { raised: !0, key: g.label, className: a.card },
      e.createElement(
        V,
        null,
        e.createElement(B, null, g.label),
        e.createElement(
          z,
          { disablePadding: !0 },
          g.values.map((y, T) => e.createElement(
            W,
            { key: T, disableGutters: !0 },
            e.createElement(b, { value: y, InputProps: {
              endAdornment: e.createElement(
                F,
                { position: "end" },
                e.createElement(
                  $,
                  { onClick: () => {
                    const S = g.values.filter((P, C) => C !== T);
                    p === 0 ? n(S) : u(S);
                  } },
                  e.createElement(K, null)
                )
              )
            } })
          )),
          e.createElement(
            W,
            { disableGutters: !0 },
            e.createElement(b, { value: p === 0 ? r : s, placeholder: "add new", onChange: (y) => {
              p === 0 ? c(y.target.value) : l(y.target.value);
            }, InputProps: {
              endAdornment: e.createElement(
                F,
                { position: "end" },
                e.createElement(
                  $,
                  { onClick: () => {
                    p === 0 ? (n([...o, r]), c("")) : (u([...i, s]), l(""));
                  }, disabled: p === 0 ? r === "" : s === "", "data-testid": "stringArrayAdd-Feat" },
                  e.createElement(Y, null)
                )
              )
            } })
          )
        )
      )
    ))
  );
}), pe = v()((t) => ({
  spacing: {
    marginBottom: t.spacing(3)
  }
})), de = k(({ model: t }) => {
  const { classes: a } = pe(), { pluginManager: r } = d.getEnv(t), { trackType: c } = t, s = r.getTrackElements();
  return e.createElement(b, { className: a.spacing, value: c, variant: "outlined", label: "Track type", helperText: "Select track type", select: !0, fullWidth: !0, onChange: (l) => t.setTrackType(l.target.value), SelectProps: {
    // @ts-expect-error
    SelectDisplayProps: { "data-testid": "trackTypeSelect" }
  } }, s.map(({ name: l, displayName: o }) => e.createElement(D, { key: l, value: l }, o)));
}), me = v()((t) => ({
  spacing: {
    marginBottom: t.spacing(3)
  }
}));
function ge(t) {
  const a = {};
  return t.forEach((r) => {
    var c;
    const s = ((c = r.adapterMetadata) === null || c === void 0 ? void 0 : c.category) || "Default";
    a[s] || (a[s] = []), a[s].push(r);
  }), a;
}
const J = k(({ model: t }) => {
  const { classes: a } = me(), { trackAdapter: r } = t, { pluginManager: c } = d.getEnv(t);
  return e.createElement(b, { className: a.spacing, value: (r == null ? void 0 : r.type) !== "UNKNOWN" ? r == null ? void 0 : r.type : "", label: "Adapter type", variant: "outlined", helperText: "Select an adapter type", select: !0, fullWidth: !0, onChange: (s) => t.setAdapterHint(s.target.value), SelectProps: {
    // @ts-expect-error
    SelectDisplayProps: { "data-testid": "adapterTypeSelect" }
  } }, Object.entries(ge(c.getAdapterElements().filter((s) => {
    var l;
    return !(!((l = s.adapterMetadata) === null || l === void 0) && l.hiddenFromGUI);
  }))).map(([s, l]) => [
    e.createElement(Z, null, s),
    l.map((o) => e.createElement(D, { key: o.name, value: o.name }, o.displayName))
  ]));
}), M = v()((t) => ({
  spacing: {
    marginBottom: t.spacing(3)
  }
}));
function fe({ trackAdapter: t, trackType: a }) {
  const { classes: r } = M(), { type: c, subadapter: s } = t;
  return c === "SNPCoverageAdapter" ? e.createElement(
    h,
    { className: r.spacing },
    "Selected ",
    e.createElement("code", null, a),
    ". Using adapter ",
    e.createElement("code", null, c),
    " with subadapter ",
    e.createElement("code", null, s == null ? void 0 : s.type),
    ". Please enter a track name and, if necessary, update the track type."
  ) : e.createElement(
    h,
    { className: r.spacing },
    "Using adapter ",
    e.createElement("code", null, c),
    " and guessing track type",
    " ",
    e.createElement("code", null, a),
    ". Please enter a track name and, if necessary, update the track type."
  );
}
function ke({ model: t }) {
  const { classes: a } = M();
  return e.createElement(
    e.Fragment,
    null,
    e.createElement(
      h,
      { className: a.spacing },
      "JBrowse was not able to guess the adapter type for this data, but it may be in the list below. If not, you can",
      " ",
      e.createElement(w, { href: "https://github.com/GMOD/jbrowse-components/releases", target: "_blank", rel: "noopener noreferrer" }, "check for new releases"),
      " ",
      "of JBrowse to see if they support this data type or",
      " ",
      e.createElement(w, { href: "https://github.com/GMOD/jbrowse-components/issues/new", target: "_blank", rel: "noopener noreferrer" }, "file an issue"),
      " ",
      "and add a feature request for this data type."
    ),
    e.createElement(J, { model: t })
  );
}
const Ee = k(function({ model: a }) {
  const { classes: r } = M(), [c, s] = m(!0), l = d.getSession(a), { trackName: o, trackAdapter: n, trackType: i, warningMessage: u, adapterHint: E } = a;
  if (L(() => {
    E === "" && n && a.setAdapterHint(n.type);
  }, [E, n, n == null ? void 0 : n.type, a]), a.unsupported)
    return e.createElement(
      h,
      { className: r.spacing },
      "This version of JBrowse cannot display data of this type. It is possible, however, that there is a newer version that can display them. You can",
      " ",
      e.createElement(w, { href: "https://github.com/GMOD/jbrowse-components/releases", target: "_blank", rel: "noopener noreferrer" }, "check for new releases"),
      " ",
      "of JBrowse or",
      " ",
      e.createElement(w, { href: "https://github.com/GMOD/jbrowse-components/issues/new", target: "_blank", rel: "noopener noreferrer" }, "file an issue"),
      " ",
      "and add a feature request for this data type."
    );
  if ((n == null ? void 0 : n.type) === Q.UNKNOWN)
    return e.createElement(ke, { model: a });
  if (!(n != null && n.type))
    return e.createElement(h, null, "Could not recognize this data type.");
  const g = d.supportedIndexingAdapters(n == null ? void 0 : n.type);
  return e.createElement(
    "div",
    null,
    n ? e.createElement(fe, { trackAdapter: n, trackType: i }) : null,
    u ? e.createElement(h, { style: { color: "orange" } }, u) : null,
    e.createElement(b, { className: r.spacing, label: "trackName", helperText: "A name for this track", fullWidth: !0, value: o, onChange: (p) => a.setTrackName(p.target.value), inputProps: { "data-testid": "trackNameInput" } }),
    e.createElement(J, { model: a }),
    e.createElement(de, { model: a }),
    e.createElement(N.AssemblySelector, { session: l, helperText: "Select assembly to add track to", selected: a.assembly, onChange: (p) => a.setAssembly(p), TextFieldProps: {
      fullWidth: !0,
      SelectProps: {
        // @ts-expect-error
        SelectDisplayProps: { "data-testid": "assemblyNameSelect" }
      }
    } }),
    d.isElectron && g && e.createElement(
      U,
      null,
      e.createElement(X, { label: "Index track for text searching?", control: e.createElement(ee, { checked: c, onChange: (p) => {
        s(p.target.checked), a.setTextIndexTrack(p.target.checked);
      } }) })
    ),
    d.isElectron && c && g ? e.createElement(ue, { model: a }) : null
  );
}), ye = v()((t) => ({
  paper: {
    padding: t.spacing(2)
  },
  spacer: {
    height: t.spacing(8)
  }
}));
function he({ model: t }) {
  const { classes: a } = ye(), r = j(t);
  return e.createElement(
    _,
    { className: a.paper },
    e.createElement(N.FileSelector, { name: "Main file", description: "", location: t.trackData, setLocation: t.setTrackData, setName: t.setTrackName, rootModel: r }),
    e.createElement("div", { className: a.spacer }),
    e.createElement(N.FileSelector, { name: "Index file", description: "(Optional) The URL of the index file is automatically inferred from the URL of the main file if it is not supplied.", location: t.indexTrackData, setLocation: t.setIndexTrackData, setName: t.setTrackName, rootModel: r })
  );
}
const ve = k(he), be = v()((t) => ({
  root: {
    marginTop: t.spacing(1)
  },
  stepper: {
    backgroundColor: t.palette.background.default
  },
  button: {
    marginRight: t.spacing(1)
  },
  actionsContainer: {
    marginTop: t.spacing(10),
    marginBottom: t.spacing(2)
  },
  alertContainer: {
    padding: `${t.spacing(2)}px 0px ${t.spacing(2)}px 0px`
  }
})), A = ["Enter track data", "Confirm track type"];
function Te({ model: t }) {
  const [a, r] = m(0), { classes: c } = be(), { jobsManager: s } = j(t), l = d.getSession(t), { assembly: o, trackAdapter: n, trackData: i, trackName: u, trackType: E, textIndexTrack: g, textIndexingConf: p } = t, [y, T] = m();
  function S(f) {
    switch (f) {
      case 0:
        return e.createElement(ve, { model: t });
      case 1:
        return e.createElement(Ee, { model: t });
      default:
        return e.createElement(h, null, "Unknown step");
    }
  }
  async function P() {
    if (a !== A.length - 1) {
      r(a + 1);
      return;
    }
    const f = [
      `${u.toLowerCase().replace(/ /g, "_")}-${Date.now()}`,
      `${l.adminMode ? "" : "-sessionTrack"}`
    ].join(""), x = l.assemblyManager.get(o);
    if (x && n && n.type !== "UNKNOWN") {
      if (l.addTrackConf({
        trackId: f,
        type: E,
        name: u,
        assemblyNames: [o],
        adapter: {
          ...n,
          sequenceAdapter: ce.getConf(x, ["sequence", "adapter"])
        }
      }), t.view) {
        if (t.view.showTrack(f), d.isElectron && g && d.supportedIndexingAdapters(n.type)) {
          const q = p || {
            attributes: ["Name", "ID"],
            exclude: ["CDS", "exon"]
          }, O = u + "-index", R = {
            indexingParams: {
              ...q,
              assemblies: [o],
              tracks: [f],
              indexType: "perTrack",
              name: O,
              timestamp: (/* @__PURE__ */ new Date()).toISOString()
            },
            name: O,
            cancelCallback: () => s.abortJob()
          };
          s.queueJob(R);
        }
      } else
        l.notify("Open a new view, or use the track selector in an existing view, to view this track", "info");
      t.clearData(), l.hideWidget(t);
    } else
      T(`Failed to add track.
The configuration of this file is not currently supported.`);
  }
  function C() {
    T(void 0), r(a - 1);
  }
  function G() {
    switch (a) {
      case 0:
        return !i;
      case 1:
        return !(u && E && (n != null && n.type) && o);
      default:
        return !0;
    }
  }
  return e.createElement(
    "div",
    { className: c.root },
    e.createElement(te, { className: c.stepper, activeStep: a, orientation: "vertical" }, A.map((f, x) => e.createElement(
      ae,
      { key: f },
      e.createElement(ne, null, f),
      e.createElement(
        re,
        null,
        S(x),
        e.createElement(
          "div",
          { className: c.actionsContainer },
          e.createElement(I, { disabled: a === 0, onClick: C, className: c.button }, "Back"),
          e.createElement(I, { disabled: G(), variant: "contained", color: "primary", onClick: P, className: c.button, "data-testid": "addTrackNextButton" }, a === A.length - 1 ? "Add" : "Next")
        ),
        y ? e.createElement(
          "div",
          { className: c.alertContainer },
          e.createElement(se, { severity: "error" }, y)
        ) : null
      )
    )))
  );
}
const Se = k(Te), xe = v()({
  textbox: {
    width: "100%"
  },
  submit: {
    marginTop: 25,
    marginBottom: 100,
    display: "block"
  }
});
function we({ model: t }) {
  const { classes: a } = xe(), [r, c] = m(""), [s, l] = m();
  return e.createElement(
    "div",
    null,
    s ? e.createElement(N.ErrorMessage, { error: s }) : null,
    e.createElement(b, { multiline: !0, rows: 10, value: r, onChange: (o) => c(o.target.value), placeholder: "Paste track config or array of track configs in JSON format", variant: "outlined", className: a.textbox }),
    e.createElement(I, { variant: "contained", className: a.submit, onClick: () => {
      try {
        l(void 0);
        const o = d.getSession(t), n = JSON.parse(r), i = Array.isArray(n) ? n : [n];
        i.forEach((u) => o.addTrackConf(u)), i.forEach((u) => u.trackId), t.clearData(), o.hideWidget(t);
      } catch (o) {
        l(o);
      }
    } }, "Submit")
  );
}
const Ne = k(we);
function Ce({ model: t }) {
  const [a, r] = m("Default add track workflow"), { pluginManager: c } = d.getEnv(t), s = c.getAddTrackWorkflowElements(), l = {
    "Default add track workflow": Se,
    "Add track JSON": Ne,
    ...Object.fromEntries(s.map((i) => [i.name, i.ReactComponent]))
  }, o = l[a] ? a : "Default add track workflow", n = l[o];
  return e.createElement(
    e.Fragment,
    null,
    e.createElement(
      U,
      null,
      e.createElement(le, { value: o, onChange: (i) => r(i.target.value) }, Object.keys(l).map((i) => e.createElement(D, { key: i, value: i }, i))),
      e.createElement(oe, null, "Type of add track workflow")
    ),
    e.createElement("br", null),
    e.createElement(n, { model: t })
  );
}
const Me = k(Ce);
export {
  Me as default
};
