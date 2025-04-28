import f, { Component as je, createElement as q, useState as w, PureComponent as ke, createRef as Ve, useRef as Ke, useMemo as he, useCallback as Ge, useEffect as Je, lazy as P, Suspense as xe } from "react";
import { aG as me, q as $, o as R, t as z, aH as Qe, p as Xe, aI as Ze, M as pe, aJ as V, aK as ae, aL as W, aM as Ye, aN as Me, aO as et, r as Oe, i as Ne, j as We, T as tt, aP as nt, aQ as it, I as B, aR as K, l as $e, a9 as rt, s as at, C as st, u as ot, aS as lt, aT as ct, n as dt, a5 as ut, B as ft } from "./index-8ee3baa4.js";
import { d as ht } from "./Clear-70a3606f.js";
var mt = typeof globalThis < "u" ? globalThis : typeof self < "u" ? self : typeof window < "u" ? window : typeof me < "u" ? me : {}, Re = {};
let k;
typeof window < "u" ? k = window : typeof self < "u" ? k = self : k = mt;
let te = null, ne = null;
const ge = 20, Q = k.clearTimeout, ve = k.setTimeout, X = k.cancelAnimationFrame || k.mozCancelAnimationFrame || k.webkitCancelAnimationFrame, _e = k.requestAnimationFrame || k.mozRequestAnimationFrame || k.webkitRequestAnimationFrame;
X == null || _e == null ? (te = Q, ne = function(e) {
  return ve(e, ge);
}) : (te = function([e, t]) {
  X(e), Q(t);
}, ne = function(e) {
  const t = _e(function() {
    Q(i), e();
  }), i = ve(function() {
    X(t), e();
  }, ge);
  return [
    t,
    i
  ];
});
function pt(r) {
  let e, t, i, a, s, o, l;
  const h = typeof document < "u" && document.attachEvent;
  if (!h) {
    o = function(d) {
      const m = d.__resizeTriggers__, g = m.firstElementChild, T = m.lastElementChild, b = g.firstElementChild;
      T.scrollLeft = T.scrollWidth, T.scrollTop = T.scrollHeight, b.style.width = g.offsetWidth + 1 + "px", b.style.height = g.offsetHeight + 1 + "px", g.scrollLeft = g.scrollWidth, g.scrollTop = g.scrollHeight;
    }, s = function(d) {
      return d.offsetWidth !== d.__resizeLast__.width || d.offsetHeight !== d.__resizeLast__.height;
    }, l = function(d) {
      if (d.target.className && typeof d.target.className.indexOf == "function" && d.target.className.indexOf("contract-trigger") < 0 && d.target.className.indexOf("expand-trigger") < 0)
        return;
      const m = this;
      o(this), this.__resizeRAF__ && te(this.__resizeRAF__), this.__resizeRAF__ = ne(function() {
        s(m) && (m.__resizeLast__.width = m.offsetWidth, m.__resizeLast__.height = m.offsetHeight, m.__resizeListeners__.forEach(function(b) {
          b.call(m, d);
        }));
      });
    };
    let u = !1, v = "";
    i = "animationstart";
    const S = "Webkit Moz O ms".split(" ");
    let n = "webkitAnimationStart animationstart oAnimationStart MSAnimationStart".split(" "), c = "";
    {
      const d = document.createElement("fakeelement");
      if (d.style.animationName !== void 0 && (u = !0), u === !1) {
        for (let m = 0; m < S.length; m++)
          if (d.style[S[m] + "AnimationName"] !== void 0) {
            c = S[m], v = "-" + c.toLowerCase() + "-", i = n[m], u = !0;
            break;
          }
      }
    }
    t = "resizeanim", e = "@" + v + "keyframes " + t + " { from { opacity: 0; } to { opacity: 0; } } ", a = v + "animation: 1ms " + t + "; ";
  }
  const p = function(u) {
    if (!u.getElementById("detectElementResize")) {
      const v = (e || "") + ".resize-triggers { " + (a || "") + 'visibility: hidden; opacity: 0; } .resize-triggers, .resize-triggers > div, .contract-trigger:before { content: " "; display: block; position: absolute; top: 0; left: 0; height: 100%; width: 100%; overflow: hidden; z-index: -1; } .resize-triggers > div { background: #eee; overflow: auto; } .contract-trigger:before { width: 200%; height: 200%; }', S = u.head || u.getElementsByTagName("head")[0], n = u.createElement("style");
      n.id = "detectElementResize", n.type = "text/css", r != null && n.setAttribute("nonce", r), n.styleSheet ? n.styleSheet.cssText = v : n.appendChild(u.createTextNode(v)), S.appendChild(n);
    }
  };
  return {
    addResizeListener: function(u, v) {
      if (h)
        u.attachEvent("onresize", v);
      else {
        if (!u.__resizeTriggers__) {
          const S = u.ownerDocument, n = k.getComputedStyle(u);
          n && n.position === "static" && (u.style.position = "relative"), p(S), u.__resizeLast__ = {}, u.__resizeListeners__ = [], (u.__resizeTriggers__ = S.createElement("div")).className = "resize-triggers";
          const c = S.createElement("div");
          c.className = "expand-trigger", c.appendChild(S.createElement("div"));
          const d = S.createElement("div");
          d.className = "contract-trigger", u.__resizeTriggers__.appendChild(c), u.__resizeTriggers__.appendChild(d), u.appendChild(u.__resizeTriggers__), o(u), u.addEventListener("scroll", l, !0), i && (u.__resizeTriggers__.__animationListener__ = function(g) {
            g.animationName === t && o(u);
          }, u.__resizeTriggers__.addEventListener(i, u.__resizeTriggers__.__animationListener__));
        }
        u.__resizeListeners__.push(v);
      }
    },
    removeResizeListener: function(u, v) {
      if (h)
        u.detachEvent("onresize", v);
      else if (u.__resizeListeners__.splice(u.__resizeListeners__.indexOf(v), 1), !u.__resizeListeners__.length) {
        u.removeEventListener("scroll", l, !0), u.__resizeTriggers__.__animationListener__ && (u.__resizeTriggers__.removeEventListener(i, u.__resizeTriggers__.__animationListener__), u.__resizeTriggers__.__animationListener__ = null);
        try {
          u.__resizeTriggers__ = !u.removeChild(u.__resizeTriggers__);
        } catch {
        }
      }
    }
  };
}
Re = {
  createDetectElementResize: pt
};
function N(r, e, t) {
  return e = gt(e), e in r ? Object.defineProperty(r, e, {
    value: t,
    enumerable: !0,
    configurable: !0,
    writable: !0
  }) : r[e] = t, r;
}
function gt(r) {
  var e = vt(r, "string");
  return typeof e == "symbol" ? e : String(e);
}
function vt(r, e) {
  if (typeof r != "object" || r === null)
    return r;
  var t = r[Symbol.toPrimitive];
  if (t !== void 0) {
    var i = t.call(r, e || "default");
    if (typeof i != "object")
      return i;
    throw new TypeError("@@toPrimitive must return a primitive value.");
  }
  return (e === "string" ? String : Number)(r);
}
class De extends je {
  constructor(...e) {
    super(...e), N(this, "state", {
      height: this.props.defaultHeight || 0,
      scaledHeight: this.props.defaultHeight || 0,
      scaledWidth: this.props.defaultWidth || 0,
      width: this.props.defaultWidth || 0
    }), N(this, "_autoSizer", null), N(this, "_detectElementResize", null), N(this, "_parentNode", null), N(this, "_resizeObserver", null), N(this, "_onResize", () => {
      const { disableHeight: t, disableWidth: i, onResize: a } = this.props;
      if (this._parentNode) {
        var s, o, l, h;
        const p = window.getComputedStyle(this._parentNode) || {}, _ = parseFloat((s = p.paddingLeft) !== null && s !== void 0 ? s : "0"), y = parseFloat((o = p.paddingRight) !== null && o !== void 0 ? o : "0"), u = parseFloat((l = p.paddingTop) !== null && l !== void 0 ? l : "0"), v = parseFloat((h = p.paddingBottom) !== null && h !== void 0 ? h : "0"), S = this._parentNode.getBoundingClientRect(), n = S.height - u - v, c = S.width - _ - y, d = this._parentNode.offsetHeight - u - v, m = this._parentNode.offsetWidth - _ - y;
        (!t && (this.state.height !== d || this.state.scaledHeight !== n) || !i && (this.state.width !== m || this.state.scaledWidth !== c)) && (this.setState({
          height: d,
          width: m,
          scaledHeight: n,
          scaledWidth: c
        }), typeof a == "function" && a({
          height: d,
          scaledHeight: n,
          scaledWidth: c,
          width: m
        }));
      }
    }), N(this, "_setRef", (t) => {
      this._autoSizer = t;
    });
  }
  componentDidMount() {
    const { nonce: e } = this.props;
    this._autoSizer && this._autoSizer.parentNode && this._autoSizer.parentNode.ownerDocument && this._autoSizer.parentNode.ownerDocument.defaultView && this._autoSizer.parentNode instanceof this._autoSizer.parentNode.ownerDocument.defaultView.HTMLElement && (this._parentNode = this._autoSizer.parentNode, this._parentNode != null && (typeof ResizeObserver < "u" ? (this._resizeObserver = new ResizeObserver(() => {
      setTimeout(this._onResize, 0);
    }), this._resizeObserver.observe(this._parentNode)) : (this._detectElementResize = (0, Re.createDetectElementResize)(e), this._detectElementResize.addResizeListener(this._parentNode, this._onResize)), this._onResize()));
  }
  componentWillUnmount() {
    this._parentNode && (this._detectElementResize && this._detectElementResize.removeResizeListener(this._parentNode, this._onResize), this._resizeObserver && (this._resizeObserver.observe(this._parentNode), this._resizeObserver.disconnect()));
  }
  render() {
    const { children: e, defaultHeight: t, defaultWidth: i, disableHeight: a, disableWidth: s, nonce: o, onResize: l, style: h, tagName: p = "div", ..._ } = this.props, { height: y, scaledHeight: u, scaledWidth: v, width: S } = this.state, n = {
      overflow: "visible"
    }, c = {};
    let d = !1;
    return a || (y === 0 && (d = !0), n.height = 0, c.height = y, c.scaledHeight = u), s || (S === 0 && (d = !0), n.width = 0, c.width = S, c.scaledWidth = v), q(p, {
      ref: this._setRef,
      style: {
        ...n,
        ...h
      },
      ..._
    }, !d && e(c));
  }
}
N(De, "defaultProps", {
  onResize: () => {
  },
  disableHeight: !1,
  disableWidth: !1,
  style: {}
});
var _t = De;
const bt = $()((r) => ({
  fab: {
    position: "absolute",
    bottom: r.spacing(6),
    right: r.spacing(6)
  }
})), St = R(function({ model: e }) {
  const { classes: t } = bt(), i = z.getSession(e), [a, s] = w(null);
  function o() {
    s(null);
  }
  const l = z.isSessionModelWithConnections(i), h = z.isSessionWithAddTracks(i);
  return f.createElement(f.Fragment, null, h || l ? f.createElement(
    f.Fragment,
    null,
    f.createElement(
      Qe,
      { color: "secondary", className: t.fab, onClick: (p) => s(p.currentTarget) },
      f.createElement(Xe, null)
    ),
    f.createElement(
      Ze,
      { anchorEl: a, open: !!a, onClose: () => s(null) },
      l ? f.createElement(pe, { onClick: () => {
        o(), z.isSessionModelWithWidgets(i) && i.showWidget(i.addWidget("AddConnectionWidget", "addConnectionWidget"));
      } }, "Add connection") : null,
      h ? f.createElement(pe, { onClick: () => {
        o(), z.isSessionModelWithWidgets(i) && i.showWidget(i.addWidget("AddTrackWidget", "addTrackWidget", {
          view: e.view.id
        }));
      } }, "Add track") : null
    )
  ) : null);
});
var be = Number.isNaN || function(e) {
  return typeof e == "number" && e !== e;
};
function Tt(r, e) {
  return !!(r === e || be(r) && be(e));
}
function yt(r, e) {
  if (r.length !== e.length)
    return !1;
  for (var t = 0; t < r.length; t++)
    if (!Tt(r[t], e[t]))
      return !1;
  return !0;
}
function Z(r, e) {
  e === void 0 && (e = yt);
  var t, i = [], a, s = !1;
  function o() {
    for (var l = [], h = 0; h < arguments.length; h++)
      l[h] = arguments[h];
    return s && t === this && e(l, i) || (a = r.apply(this, l), s = !0, t = this, i = l), a;
  }
  return o;
}
var It = typeof performance == "object" && typeof performance.now == "function", Se = It ? function() {
  return performance.now();
} : function() {
  return Date.now();
};
function Te(r) {
  cancelAnimationFrame(r.id);
}
function Et(r, e) {
  var t = Se();
  function i() {
    Se() - t >= e ? r.call(null) : a.id = requestAnimationFrame(i);
  }
  var a = {
    id: requestAnimationFrame(i)
  };
  return a;
}
var Y = -1;
function ye(r) {
  if (r === void 0 && (r = !1), Y === -1 || r) {
    var e = document.createElement("div"), t = e.style;
    t.width = "50px", t.height = "50px", t.overflow = "scroll", document.body.appendChild(e), Y = e.offsetWidth - e.clientWidth, document.body.removeChild(e);
  }
  return Y;
}
var L = null;
function Ie(r) {
  if (r === void 0 && (r = !1), L === null || r) {
    var e = document.createElement("div"), t = e.style;
    t.width = "50px", t.height = "50px", t.overflow = "scroll", t.direction = "rtl";
    var i = document.createElement("div"), a = i.style;
    return a.width = "100px", a.height = "100px", e.appendChild(i), document.body.appendChild(e), e.scrollLeft > 0 ? L = "positive-descending" : (e.scrollLeft = 1, e.scrollLeft === 0 ? L = "negative" : L = "positive-ascending"), document.body.removeChild(e), L;
  }
  return L;
}
V.env.NODE_ENV;
var zt = 150, Ct = function(e, t) {
  return e;
}, U = null, j = null;
V.env.NODE_ENV !== "production" && typeof window < "u" && typeof window.WeakSet < "u" && (U = /* @__PURE__ */ new WeakSet(), j = /* @__PURE__ */ new WeakSet());
function wt(r) {
  var e, t = r.getItemOffset, i = r.getEstimatedTotalSize, a = r.getItemSize, s = r.getOffsetForIndexAndAlignment, o = r.getStartIndexForOffset, l = r.getStopIndexForStartIndex, h = r.initInstanceProps, p = r.shouldResetStyleCacheOnItemSizeChange, _ = r.validateProps;
  return e = /* @__PURE__ */ function(y) {
    ae(u, y);
    function u(S) {
      var n;
      return n = y.call(this, S) || this, n._instanceProps = h(n.props, W(n)), n._outerRef = void 0, n._resetIsScrollingTimeoutId = null, n.state = {
        instance: W(n),
        isScrolling: !1,
        scrollDirection: "forward",
        scrollOffset: typeof n.props.initialScrollOffset == "number" ? n.props.initialScrollOffset : 0,
        scrollUpdateWasRequested: !1
      }, n._callOnItemsRendered = void 0, n._callOnItemsRendered = Z(function(c, d, m, g) {
        return n.props.onItemsRendered({
          overscanStartIndex: c,
          overscanStopIndex: d,
          visibleStartIndex: m,
          visibleStopIndex: g
        });
      }), n._callOnScroll = void 0, n._callOnScroll = Z(function(c, d, m) {
        return n.props.onScroll({
          scrollDirection: c,
          scrollOffset: d,
          scrollUpdateWasRequested: m
        });
      }), n._getItemStyle = void 0, n._getItemStyle = function(c) {
        var d = n.props, m = d.direction, g = d.itemSize, T = d.layout, b = n._getItemStyleCache(p && g, p && T, p && m), I;
        if (b.hasOwnProperty(c))
          I = b[c];
        else {
          var E = t(n.props, c, n._instanceProps), x = a(n.props, c, n._instanceProps), C = m === "horizontal" || T === "horizontal", M = m === "rtl", D = C ? E : 0;
          b[c] = I = {
            position: "absolute",
            left: M ? void 0 : D,
            right: M ? D : void 0,
            top: C ? 0 : E,
            height: C ? "100%" : x,
            width: C ? x : "100%"
          };
        }
        return I;
      }, n._getItemStyleCache = void 0, n._getItemStyleCache = Z(function(c, d, m) {
        return {};
      }), n._onScrollHorizontal = function(c) {
        var d = c.currentTarget, m = d.clientWidth, g = d.scrollLeft, T = d.scrollWidth;
        n.setState(function(b) {
          if (b.scrollOffset === g)
            return null;
          var I = n.props.direction, E = g;
          if (I === "rtl")
            switch (Ie()) {
              case "negative":
                E = -g;
                break;
              case "positive-descending":
                E = T - m - g;
                break;
            }
          return E = Math.max(0, Math.min(E, T - m)), {
            isScrolling: !0,
            scrollDirection: b.scrollOffset < g ? "forward" : "backward",
            scrollOffset: E,
            scrollUpdateWasRequested: !1
          };
        }, n._resetIsScrollingDebounced);
      }, n._onScrollVertical = function(c) {
        var d = c.currentTarget, m = d.clientHeight, g = d.scrollHeight, T = d.scrollTop;
        n.setState(function(b) {
          if (b.scrollOffset === T)
            return null;
          var I = Math.max(0, Math.min(T, g - m));
          return {
            isScrolling: !0,
            scrollDirection: b.scrollOffset < I ? "forward" : "backward",
            scrollOffset: I,
            scrollUpdateWasRequested: !1
          };
        }, n._resetIsScrollingDebounced);
      }, n._outerRefSetter = function(c) {
        var d = n.props.outerRef;
        n._outerRef = c, typeof d == "function" ? d(c) : d != null && typeof d == "object" && d.hasOwnProperty("current") && (d.current = c);
      }, n._resetIsScrollingDebounced = function() {
        n._resetIsScrollingTimeoutId !== null && Te(n._resetIsScrollingTimeoutId), n._resetIsScrollingTimeoutId = Et(n._resetIsScrolling, zt);
      }, n._resetIsScrolling = function() {
        n._resetIsScrollingTimeoutId = null, n.setState({
          isScrolling: !1
        }, function() {
          n._getItemStyleCache(-1, null);
        });
      }, n;
    }
    u.getDerivedStateFromProps = function(n, c) {
      return kt(n, c), _(n), null;
    };
    var v = u.prototype;
    return v.scrollTo = function(n) {
      n = Math.max(0, n), this.setState(function(c) {
        return c.scrollOffset === n ? null : {
          scrollDirection: c.scrollOffset < n ? "forward" : "backward",
          scrollOffset: n,
          scrollUpdateWasRequested: !0
        };
      }, this._resetIsScrollingDebounced);
    }, v.scrollToItem = function(n, c) {
      c === void 0 && (c = "auto");
      var d = this.props, m = d.itemCount, g = d.layout, T = this.state.scrollOffset;
      n = Math.max(0, Math.min(n, m - 1));
      var b = 0;
      if (this._outerRef) {
        var I = this._outerRef;
        g === "vertical" ? b = I.scrollWidth > I.clientWidth ? ye() : 0 : b = I.scrollHeight > I.clientHeight ? ye() : 0;
      }
      this.scrollTo(s(this.props, n, c, T, this._instanceProps, b));
    }, v.componentDidMount = function() {
      var n = this.props, c = n.direction, d = n.initialScrollOffset, m = n.layout;
      if (typeof d == "number" && this._outerRef != null) {
        var g = this._outerRef;
        c === "horizontal" || m === "horizontal" ? g.scrollLeft = d : g.scrollTop = d;
      }
      this._callPropsCallbacks();
    }, v.componentDidUpdate = function() {
      var n = this.props, c = n.direction, d = n.layout, m = this.state, g = m.scrollOffset, T = m.scrollUpdateWasRequested;
      if (T && this._outerRef != null) {
        var b = this._outerRef;
        if (c === "horizontal" || d === "horizontal")
          if (c === "rtl")
            switch (Ie()) {
              case "negative":
                b.scrollLeft = -g;
                break;
              case "positive-ascending":
                b.scrollLeft = g;
                break;
              default:
                var I = b.clientWidth, E = b.scrollWidth;
                b.scrollLeft = E - I - g;
                break;
            }
          else
            b.scrollLeft = g;
        else
          b.scrollTop = g;
      }
      this._callPropsCallbacks();
    }, v.componentWillUnmount = function() {
      this._resetIsScrollingTimeoutId !== null && Te(this._resetIsScrollingTimeoutId);
    }, v.render = function() {
      var n = this.props, c = n.children, d = n.className, m = n.direction, g = n.height, T = n.innerRef, b = n.innerElementType, I = n.innerTagName, E = n.itemCount, x = n.itemData, C = n.itemKey, M = C === void 0 ? Ct : C, D = n.layout, G = n.outerElementType, A = n.outerTagName, O = n.style, He = n.useIsScrolling, Pe = n.width, ce = this.state.isScrolling, J = m === "horizontal" || D === "horizontal", Be = J ? this._onScrollHorizontal : this._onScrollVertical, de = this._getRangeToRender(), qe = de[0], Ue = de[1], ue = [];
      if (E > 0)
        for (var H = qe; H <= Ue; H++)
          ue.push(q(c, {
            data: x,
            key: M(H, x),
            index: H,
            isScrolling: He ? ce : void 0,
            style: this._getItemStyle(H)
          }));
      var fe = i(this.props, this._instanceProps);
      return q(G || A || "div", {
        className: d,
        onScroll: Be,
        ref: this._outerRefSetter,
        style: Ye({
          position: "relative",
          height: g,
          width: Pe,
          overflow: "auto",
          WebkitOverflowScrolling: "touch",
          willChange: "transform",
          direction: m
        }, O)
      }, q(b || I || "div", {
        children: ue,
        ref: T,
        style: {
          height: J ? "100%" : fe,
          pointerEvents: ce ? "none" : void 0,
          width: J ? fe : "100%"
        }
      }));
    }, v._callPropsCallbacks = function() {
      if (typeof this.props.onItemsRendered == "function") {
        var n = this.props.itemCount;
        if (n > 0) {
          var c = this._getRangeToRender(), d = c[0], m = c[1], g = c[2], T = c[3];
          this._callOnItemsRendered(d, m, g, T);
        }
      }
      if (typeof this.props.onScroll == "function") {
        var b = this.state, I = b.scrollDirection, E = b.scrollOffset, x = b.scrollUpdateWasRequested;
        this._callOnScroll(I, E, x);
      }
    }, v._getRangeToRender = function() {
      var n = this.props, c = n.itemCount, d = n.overscanCount, m = this.state, g = m.isScrolling, T = m.scrollDirection, b = m.scrollOffset;
      if (c === 0)
        return [0, 0, 0, 0];
      var I = o(this.props, b, this._instanceProps), E = l(this.props, I, b, this._instanceProps), x = !g || T === "backward" ? Math.max(1, d) : 1, C = !g || T === "forward" ? Math.max(1, d) : 1;
      return [Math.max(0, I - x), Math.max(0, Math.min(c - 1, E + C)), I, E];
    }, u;
  }(ke), e.defaultProps = {
    direction: "ltr",
    itemData: void 0,
    layout: "vertical",
    overscanCount: 2,
    useIsScrolling: !1
  }, e;
}
var kt = function(e, t) {
  var i = e.children, a = e.direction, s = e.height, o = e.layout, l = e.innerTagName, h = e.outerTagName, p = e.width, _ = t.instance;
  if (V.env.NODE_ENV !== "production") {
    (l != null || h != null) && j && !j.has(_) && (j.add(_), console.warn("The innerTagName and outerTagName props have been deprecated. Please use the innerElementType and outerElementType props instead."));
    var y = a === "horizontal" || o === "horizontal";
    switch (a) {
      case "horizontal":
      case "vertical":
        U && !U.has(_) && (U.add(_), console.warn('The direction prop should be either "ltr" (default) or "rtl". Please use the layout prop to specify "vertical" (default) or "horizontal" orientation.'));
        break;
      case "ltr":
      case "rtl":
        break;
      default:
        throw Error('An invalid "direction" prop has been specified. Value should be either "ltr" or "rtl". ' + ('"' + a + '" was specified.'));
    }
    switch (o) {
      case "horizontal":
      case "vertical":
        break;
      default:
        throw Error('An invalid "layout" prop has been specified. Value should be either "horizontal" or "vertical". ' + ('"' + o + '" was specified.'));
    }
    if (i == null)
      throw Error('An invalid "children" prop has been specified. Value should be a React component. ' + ('"' + (i === null ? "null" : typeof i) + '" was specified.'));
    if (y && typeof p != "number")
      throw Error('An invalid "width" prop has been specified. Horizontal lists must specify a number for width. ' + ('"' + (p === null ? "null" : typeof p) + '" was specified.'));
    if (!y && typeof s != "number")
      throw Error('An invalid "height" prop has been specified. Vertical lists must specify a number for height. ' + ('"' + (s === null ? "null" : typeof s) + '" was specified.'));
  }
}, xt = 50, F = function(e, t, i) {
  var a = e, s = a.itemSize, o = i.itemMetadataMap, l = i.lastMeasuredIndex;
  if (t > l) {
    var h = 0;
    if (l >= 0) {
      var p = o[l];
      h = p.offset + p.size;
    }
    for (var _ = l + 1; _ <= t; _++) {
      var y = s(_);
      o[_] = {
        offset: h,
        size: y
      }, h += y;
    }
    i.lastMeasuredIndex = t;
  }
  return o[t];
}, Mt = function(e, t, i) {
  var a = t.itemMetadataMap, s = t.lastMeasuredIndex, o = s > 0 ? a[s].offset : 0;
  return o >= i ? Ae(e, t, s, 0, i) : Ot(e, t, Math.max(0, s), i);
}, Ae = function(e, t, i, a, s) {
  for (; a <= i; ) {
    var o = a + Math.floor((i - a) / 2), l = F(e, o, t).offset;
    if (l === s)
      return o;
    l < s ? a = o + 1 : l > s && (i = o - 1);
  }
  return a > 0 ? a - 1 : 0;
}, Ot = function(e, t, i, a) {
  for (var s = e.itemCount, o = 1; i < s && F(e, i, t).offset < a; )
    i += o, o *= 2;
  return Ae(e, t, Math.min(i, s - 1), Math.floor(i / 2), a);
}, Ee = function(e, t) {
  var i = e.itemCount, a = t.itemMetadataMap, s = t.estimatedItemSize, o = t.lastMeasuredIndex, l = 0;
  if (o >= i && (o = i - 1), o >= 0) {
    var h = a[o];
    l = h.offset + h.size;
  }
  var p = i - o - 1, _ = p * s;
  return l + _;
}, Nt = /* @__PURE__ */ wt({
  getItemOffset: function(e, t, i) {
    return F(e, t, i).offset;
  },
  getItemSize: function(e, t, i) {
    return i.itemMetadataMap[t].size;
  },
  getEstimatedTotalSize: Ee,
  getOffsetForIndexAndAlignment: function(e, t, i, a, s, o) {
    var l = e.direction, h = e.height, p = e.layout, _ = e.width, y = l === "horizontal" || p === "horizontal", u = y ? _ : h, v = F(e, t, s), S = Ee(e, s), n = Math.max(0, Math.min(S - u, v.offset)), c = Math.max(0, v.offset - u + v.size + o);
    switch (i === "smart" && (a >= c - u && a <= n + u ? i = "auto" : i = "center"), i) {
      case "start":
        return n;
      case "end":
        return c;
      case "center":
        return Math.round(c + (n - c) / 2);
      case "auto":
      default:
        return a >= c && a <= n ? a : a < c ? c : n;
    }
  },
  getStartIndexForOffset: function(e, t, i) {
    return Mt(e, i, t);
  },
  getStopIndexForStartIndex: function(e, t, i, a) {
    for (var s = e.direction, o = e.height, l = e.itemCount, h = e.layout, p = e.width, _ = s === "horizontal" || h === "horizontal", y = _ ? p : o, u = F(e, t, a), v = i + y, S = u.offset + u.size, n = t; n < l - 1 && S < v; )
      n++, S += F(e, n, a).size;
    return n;
  },
  initInstanceProps: function(e, t) {
    var i = e, a = i.estimatedItemSize, s = {
      itemMetadataMap: {},
      estimatedItemSize: a || xt,
      lastMeasuredIndex: -1
    };
    return t.resetAfterIndex = function(o, l) {
      l === void 0 && (l = !0), s.lastMeasuredIndex = Math.min(s.lastMeasuredIndex, o - 1), t._getItemStyleCache(-1), l && t.forceUpdate();
    }, s;
  },
  shouldResetStyleCacheOnItemSizeChange: !1,
  validateProps: function(e) {
    var t = e.itemSize;
    if (V.env.NODE_ENV !== "production" && typeof t != "function")
      throw Error('An invalid "itemSize" prop has been specified. Value should be a function. ' + ('"' + (t === null ? "null" : typeof t) + '" was specified.'));
  }
});
function Wt(r) {
  return function(e) {
    r.forEach(function(t) {
      typeof t == "function" ? t(e) : t != null && (t.current = e);
    });
  };
}
var ee = function() {
}, $t = function(e, t) {
  return t === void 0 && (t = null), {
    child: null,
    isShown: t ? t.public.isOpen && t.isShown : !0,
    parent: t,
    public: e,
    sibling: null,
    visited: !1
  };
}, Rt = function(e, t) {
  var i = t.getRecordData, a = i(e), s = a.data.id;
  return s;
}, Dt = function(e) {
  var t = e.index, i = e.data, a = i.component, s = i.getRecordData, o = i.treeData, l = e.style, h = e.isScrolling, p = s(t);
  return /* @__PURE__ */ f.createElement(a, Object.assign({
    isScrolling: h,
    style: l,
    treeData: o
  }, p));
}, At = function(e, t, i) {
  var a = e.createRecord, s = t.buildingTaskTimeout, o = t.placeholder, l = t.async, h = l === void 0 ? !1 : l, p = t.treeWalker, _ = (
    // eslint-disable-next-line @typescript-eslint/no-unnecessary-condition
    h && i.records !== void 0
  ), y = i.records, u = [], v = /* @__PURE__ */ new Map(), S = s ? {
    timeout: s
  } : void 0, n = /* @__PURE__ */ new WeakMap(), c = p(), d = c.next(), m = d.value, g = a(m.data, i, void 0, _ ? y.get(m.data.id) : void 0);
  v.set(g.public.data.id, g), n.set(g, m);
  var T = g, b = !0, I = g, E = o !== void 0 && // If placeholder is set to null and this is the first build, idle callback
  // won't be used. It is necessary for trees with async data which can be
  // extremely complex but the first build is quite easy. During the following
  // idle callbacks the old tree will be shown.
  !(o === null && !i.order), x = E ? function(M) {
    return M.timeRemaining() > 0;
  } : function() {
    return !0;
  }, C = function M(D) {
    for (; T !== null; ) {
      if (!x(D)) {
        requestIdleCallback(M, S);
        return;
      }
      if (T.visited)
        T.visited = !1, T = T.sibling !== null ? T.sibling : T.parent, I = T;
      else {
        var G = c.next(n.get(T)), A = G.value;
        if (A === void 0) {
          b ? b = !1 : (T.isShown && u.push(T.public.data.id), T.visited = T.child !== null, T = T.child !== null ? T.child : T.sibling !== null ? T.sibling : T.parent), I = T;
          continue;
        }
        var O = a(A.data, i, b ? void 0 : T, _ ? y.get(A.data.id) : void 0);
        v.set(O.public.data.id, O), n.set(O, A), !b && I === T ? I.child = O : I.sibling = O, I = O;
      }
    }
    E && i.setState({
      order: u,
      records: v,
      updateRequest: {}
    });
  };
  return E ? requestIdleCallback(C, S) : C(), o !== void 0 && h && i.order ? i : {
    order: u,
    records: v
  };
}, ze = 32768, Lt = 2, Ft = function(e, t) {
  var i = e.order, a = e.records, s = t.opennessState;
  if (typeof s != "object")
    return null;
  for (var o in s)
    if (a.has(o)) {
      var l = s[o], h = a.get(o), p = typeof l == "boolean" ? {
        open: l
      } : l, _ = p.open, y = p.subtreeCallback, u = y === void 0 ? ee : y, v = ee, S = ee;
      h.isShown && (_ ? function() {
        for (var c = i.indexOf(o), d = h; d !== null; ) {
          if (d.sibling !== null) {
            d = d.sibling;
            break;
          }
          d = d.parent;
        }
        var m = d === null ? i.length - 1 - c : i.indexOf(d.public.data.id) - 1 - c, g = [[c + 1, m]];
        v = function(b) {
          if (b.isShown = b.parent ? b.parent.public.isOpen && b.parent.isShown : !0, b.isShown) {
            var I = g[g.length - 1];
            I.push(b.public.data.id), I.length === ze + Lt && g.push([c + 1 + ze * g.length, 0]);
          }
        }, S = function() {
          for (var b = 0; b < g.length; b++) {
            var I;
            (I = i).splice.apply(I, g[b]);
          }
        };
      }() : h.public.isOpen && function() {
        var c = i.indexOf(o), d = 0;
        v = function(g) {
          g.isShown && (d += 1), g.isShown = g.parent ? g.parent.public.isOpen && g.parent.isShown : !0;
        }, S = function() {
          i.splice(c + 1, d);
        };
      }());
      for (var n = h; n !== null; )
        n.visited ? (n.visited = !1, n = n === h ? null : n.sibling !== null ? n.sibling : n.parent) : (n.public.isOpen = n === h ? _ : n.public.isOpen, u(n.public, h.public), n !== h && v(n), n.visited = n.child !== null, n = // Look for child in any case
        n.child !== null ? n.child : (
          // Stop looking for next element if currentRecord is root.
          n === h ? null : (
            // Otherwise, look for sibling or parent
            n.sibling !== null ? n.sibling : n.parent
          )
        ));
      S();
    }
  return {
    order: i,
    records: a,
    updateRequest: {}
  };
}, Ht = function(e) {
  return function(t, i, a) {
    return a.refresh ? At(e, t, i) : Ft(i, a);
  };
}, Le = /* @__PURE__ */ function(r) {
  ae(e, r), e.getDerivedStateFromProps = function(a, s) {
    var o = a.listRef, l = o === void 0 ? null : o, h = a.treeWalker, p = s.computeTree, _ = s.list, y = s.order, u = s.treeWalker;
    return Me({
      attachRefs: Wt([_, l])
    }, h !== u || !y ? p(a, s, {
      refresh: !0
    }) : null, {
      treeWalker: h
    });
  };
  function e(i, a) {
    var s;
    return s = r.call(this, i, a) || this, s.getRecordData = s.getRecordData.bind(W(s)), s.state = {
      list: /* @__PURE__ */ Ve(),
      recomputeTree: s.recomputeTree.bind(W(s)),
      setState: s.setState.bind(W(s))
    }, s;
  }
  var t = e.prototype;
  return t.getItemData = function() {
    var a = this.props, s = a.children, o = a.itemData;
    return {
      component: s,
      // eslint-disable-next-line @typescript-eslint/unbound-method
      getRecordData: this.getRecordData,
      treeData: o
    };
  }, t.getRecordData = function(a) {
    var s = this.state, o = s.order, l = s.records;
    return l.get(o[a]).public;
  }, t.recomputeTree = function(a) {
    var s = this;
    return new Promise(function(o) {
      s.setState(function(l) {
        return l.computeTree(s.props, l, {
          opennessState: a
        });
      }, o);
    });
  }, t.scrollTo = function(a) {
    var s;
    (s = this.state.list.current) == null || s.scrollTo(a);
  }, t.scrollToItem = function(a, s) {
    var o;
    (o = this.state.list.current) == null || o.scrollToItem(this.state.order.indexOf(a), s);
  }, e;
}(ke);
Le.defaultProps = {
  rowComponent: Dt
};
var Pt = Ht({
  createRecord: function(e, t, i, a) {
    var s = t.recomputeTree, o = t.resetAfterId, l = $t({
      data: e,
      height: a ? a.public.height : e.defaultHeight,
      isOpen: a ? a.public.isOpen : e.isOpenByDefault,
      resize: function(p, _) {
        l.public.height = p, o(l.public.data.id, _);
      },
      setOpen: function(p) {
        var _;
        return s((_ = {}, _[e.id] = p, _));
      }
    }, i);
    return l;
  }
}), Bt = /* @__PURE__ */ function(r) {
  ae(e, r);
  function e(i, a) {
    var s;
    return s = r.call(this, i, a) || this, s.getItemSize = s.getItemSize.bind(W(s)), s.state = Me({}, s.state, {
      computeTree: Pt,
      resetAfterId: s.resetAfterId.bind(W(s))
    }), s;
  }
  var t = e.prototype;
  return t.resetAfterId = function(a, s) {
    var o;
    s === void 0 && (s = !1);
    var l = this.state, h = l.list, p = l.order;
    (o = h.current) == null || o.resetAfterIndex(p.indexOf(a), s);
  }, t.recomputeTree = function(a) {
    var s = this;
    return r.prototype.recomputeTree.call(this, a).then(function() {
      var o;
      (o = s.state.list.current) == null || o.resetAfterIndex(0, !0);
    });
  }, t.render = function() {
    var a = this.props;
    a.children;
    var s = a.placeholder, o = a.itemSize, l = a.rowComponent;
    a.treeWalker;
    var h = et(a, ["children", "placeholder", "itemSize", "rowComponent", "treeWalker"]), p = this.state, _ = p.attachRefs, y = p.order;
    return s && y.length === 0 ? s : /* @__PURE__ */ f.createElement(Nt, Object.assign({}, h, {
      itemCount: y.length,
      itemData: this.getItemData(),
      itemKey: Rt,
      itemSize: o ?? this.getItemSize,
      ref: _
    }), l);
  }, t.getItemSize = function(a) {
    return this.getRecordData(a).height;
  }, e;
}(Le), se = {}, qt = Ne;
Object.defineProperty(se, "__esModule", {
  value: !0
});
var oe = se.default = void 0, Ut = qt(Oe()), jt = We, Vt = (0, Ut.default)(/* @__PURE__ */ (0, jt.jsx)("path", {
  d: "M6 10c-1.1 0-2 .9-2 2s.9 2 2 2 2-.9 2-2-.9-2-2-2zm12 0c-1.1 0-2 .9-2 2s.9 2 2 2 2-.9 2-2-.9-2-2-2zm-6 0c-1.1 0-2 .9-2 2s.9 2 2 2 2-.9 2-2-.9-2-2-2z"
}), "MoreHoriz");
oe = se.default = Vt;
function ie(r) {
  return ((r == null ? void 0 : r.children.map((e) => e.children.length ? ie(e) : e.conf)) || []).flat(1 / 0);
}
function re(r, e = /* @__PURE__ */ new Map()) {
  r.id && r.children.length && e.set(r.id, r);
  for (let t = 0; t < r.children.length; t++) {
    const i = r.children[t];
    re(i, e);
  }
  return e;
}
function Kt(r = "") {
  return r.endsWith("(Unsupported)") || r.endsWith("(Unknown)");
}
const Gt = $()((r) => ({
  contrastColor: {
    color: r.palette.tertiary.contrastText
  },
  // margin:auto 0 to center text vertically
  accordionText: {
    margin: "auto 0",
    // width 100 so you can click anywhere on the category bar
    width: "100%"
  }
}));
function Jt({ isOpen: r, setOpen: e, data: t }) {
  const { classes: i } = Gt(), [a, s] = w(null), { name: o, model: l, id: h, tree: p, toggleCollapse: _ } = t;
  return f.createElement(
    "div",
    { className: i.accordionText, onClick: () => {
      a || (_(h), e(!r));
    } },
    f.createElement(
      tt,
      null,
      r ? f.createElement(nt, null) : f.createElement(it, null),
      o,
      f.createElement(
        B,
        { onClick: (y) => {
          s(y.currentTarget), y.stopPropagation();
        }, className: i.contrastColor },
        f.createElement(oe, null)
      )
    ),
    a ? f.createElement(K, { anchorEl: a, menuItems: [
      {
        label: "Add to selection",
        onClick: () => {
          const y = re(p).get(h);
          l.addToSelection(ie(y));
        }
      },
      {
        label: "Remove from selection",
        onClick: () => {
          const y = re(p).get(h);
          l.removeFromSelection(ie(y));
        }
      }
    ], onMenuItemClick: (y, u) => {
      u(), s(null);
    }, open: !!a, onClose: () => s(null) }) : null
  );
}
const Qt = $()((r) => ({
  compactCheckbox: {
    padding: 0
  },
  checkboxLabel: {
    marginRight: 0,
    "&:hover": {
      backgroundColor: r.palette.action.selected
    }
  }
}));
function Xt({ data: r }) {
  var e, t;
  const { classes: i } = Qt(), { checked: a, conf: s, model: o, drawerPosition: l, id: h, name: p, onChange: _, selected: y } = r, [u, v] = w(), S = s && $e.readConfObject(s, ["description"]) || "";
  return f.createElement(
    f.Fragment,
    null,
    f.createElement(
      rt,
      { title: S + (y ? " (in selection)" : ""), placement: l === "left" ? "right" : "left" },
      f.createElement(at, { className: i.checkboxLabel, control: f.createElement(st, { className: i.compactCheckbox, checked: a, onChange: () => _(h), disabled: Kt(p), inputProps: {
        // @ts-expect-error
        "data-testid": `htsTrackEntry-${h}`
      } }), label: f.createElement(
        "div",
        { style: { background: y ? "#cccc" : void 0 } },
        f.createElement(ot.SanitizedHTML, { html: p })
      ) })
    ),
    f.createElement(
      B,
      { onClick: (n) => v({ target: n.currentTarget, id: h, conf: s }), style: { padding: 0 }, "data-testid": `htsTrackEntryMenu-${h}` },
      f.createElement(oe, null)
    ),
    u ? f.createElement(K, { anchorEl: u == null ? void 0 : u.target, menuItems: [
      ...((t = (e = z.getSession(o)).getTrackActionMenuItems) === null || t === void 0 ? void 0 : t.call(e, u.conf)) || [],
      {
        label: "Add to selection",
        onClick: () => o.addToSelection([u.conf])
      },
      ...y ? [
        {
          label: "Remove from selection",
          onClick: () => o.removeFromSelection([u.conf])
        }
      ] : []
    ], onMenuItemClick: (n, c) => {
      c(), v(void 0);
    }, open: !!u, onClose: () => v(void 0) }) : null
  );
}
const Zt = $()((r) => ({
  // this accordionBase element's small padding is used to give a margin to
  // accordionColor it a "margin" because the virtualized elements can't really
  // use margin in a conventional way (it doesn't affect layout)
  accordionBase: {
    display: "flex"
  },
  accordionCard: {
    padding: 3,
    cursor: "pointer",
    display: "flex"
  },
  nestingLevelMarker: {
    position: "absolute",
    borderLeft: "1.5px solid #555"
  },
  // accordionColor set's display:flex so that the child accordionText use
  // vertically centered text
  accordionColor: {
    background: r.palette.tertiary.main,
    color: r.palette.tertiary.contrastText,
    width: "100%",
    display: "flex",
    paddingLeft: 5
  }
}));
function Yt({ data: r, isOpen: e, style: t, setOpen: i }) {
  const { isLeaf: a, nestingLevel: s } = r, { classes: o } = Zt(), l = 10, h = s * l + (a ? l : 0);
  return f.createElement(
    "div",
    { style: t, className: a ? void 0 : o.accordionBase },
    new Array(s).fill(0).map((p, _) => f.createElement("div", { key: `mark-${_}`, style: { left: _ * l + 4, height: t == null ? void 0 : t.height }, className: o.nestingLevelMarker })),
    f.createElement(
      "div",
      { className: a ? void 0 : o.accordionCard, style: {
        marginLeft: h,
        whiteSpace: "nowrap",
        width: "100%"
      } },
      f.createElement("div", { className: a ? void 0 : o.accordionColor }, a ? f.createElement(Xt, { data: r }) : f.createElement(Jt, { isOpen: e, data: r, setOpen: i }))
    )
  );
}
function Ce(r, e, t, i) {
  const a = !!r.conf, s = !!i[r.id];
  return {
    data: {
      defaultHeight: a ? 22 : 40,
      isLeaf: a,
      isOpenByDefault: !0,
      nestingLevel: e,
      selected: s,
      ...r,
      ...t
    },
    nestingLevel: e,
    node: r
  };
}
const we = R(function({ height: e, tree: t, model: i }) {
  const { filterText: a, selection: s, view: o } = i, l = Ke(null), h = z.getSession(i), { drawerPosition: p } = h, _ = he(() => Object.fromEntries(s.map((v) => [v.trackId, v])), [s]), y = he(() => ({
    onChange: (v) => o.toggleTrack(v),
    toggleCollapse: (v) => i.toggleCategory(v),
    tree: t,
    model: i,
    drawerPosition: p
  }), [o, i, p, t]), u = Ge(function* () {
    for (let S = 0; S < t.children.length; S++) {
      const n = t.children[S];
      yield Ce(n, 0, y, _);
    }
    for (; ; ) {
      const S = yield;
      for (let n = 0; n < S.node.children.length; n++) {
        const c = S.node.children[n];
        yield Ce(c, S.nestingLevel + 1, y, _);
      }
    }
  }, [t, y, _]);
  return Je(() => {
    l.current.recomputeTree({
      refreshNodes: !0,
      useDefaultHeight: !0
    });
  }, [t, a]), f.createElement(
    f.Fragment,
    null,
    f.createElement(Bt, { ref: l, treeWalker: u, height: e }, Yt)
  );
}), en = P(() => import("./CloseConnectionDialog-efd48f75.js")), tn = P(() => import("./DeleteConnectionDialog-00ad326d.js")), nn = P(() => import("./ManageConnectionsDialog-23d07bdc.js")), rn = P(() => import("./ToggleConnectionsDialog-c75042f9.js")), an = $()((r) => ({
  menuIcon: {
    marginRight: r.spacing(1),
    marginBottom: 0
  }
})), sn = R(function({ model: e, setAssemblyIdx: t }) {
  const i = z.getSession(e), [a, s] = w(), [o, l] = w(), [h, p] = w(), [_, y] = w(!1), [u, v] = w(!1), { classes: S } = an(), { assemblyNames: n } = e;
  function c(m, g) {
    const T = $e.readConfObject(m, "name"), b = i.prepareToBreakConnection(m);
    if (b) {
      const [I, E] = b;
      Object.keys(E).length > 0 ? l({
        connectionConf: m,
        safelyBreakConnection: I,
        dereferenceTypeCount: E,
        name: T
      }) : I();
    }
    g && p({ name: T, connectionConf: m });
  }
  const d = [
    {
      label: "Turn on/off connections...",
      onClick: () => y(!0)
    }
  ];
  return z.isSessionModelWithConnections(i) && (d.unshift({
    label: "Add connection...",
    onClick: () => {
      z.isSessionModelWithWidgets(i) && i.showWidget(i.addWidget("AddConnectionWidget", "addConnectionWidget"));
    }
  }), d.push({
    label: "Delete connections...",
    onClick: () => v(!0)
  })), f.createElement(
    f.Fragment,
    null,
    f.createElement(
      B,
      { className: S.menuIcon, onClick: (m) => s(m.currentTarget) },
      f.createElement(lt, null)
    ),
    f.createElement(K, { anchorEl: a, open: !!a, onMenuItemClick: (m, g) => {
      g(), s(void 0);
    }, onClose: () => s(void 0), menuItems: [
      {
        label: "Add track...",
        onClick: () => {
          z.isSessionModelWithWidgets(i) && i.showWidget(i.addWidget("AddTrackWidget", "addTrackWidget", {
            view: e.view.id
          }));
        }
      },
      ...i.makeConnection ? d : [],
      ...n.length > 1 ? [
        {
          label: "Select assembly...",
          subMenu: n.map((m, g) => ({
            label: m,
            onClick: () => t(g)
          }))
        }
      ] : []
    ] }),
    f.createElement(
      xe,
      { fallback: f.createElement("div", null) },
      o ? f.createElement(en, { modalInfo: o, setModalInfo: l }) : null,
      h ? f.createElement(tn, { handleClose: () => p(void 0), deleteDialogDetails: h, session: i }) : null,
      u ? f.createElement(nn, { handleClose: () => v(!1), breakConnection: c, session: i }) : null,
      _ ? f.createElement(rn, { handleClose: () => y(!1), session: i, breakConnection: c }) : null
    )
  );
});
var le = {}, on = Ne;
Object.defineProperty(le, "__esModule", {
  value: !0
});
var Fe = le.default = void 0, ln = on(Oe()), cn = We, dn = (0, ln.default)(/* @__PURE__ */ (0, cn.jsx)("path", {
  d: "M7 18c-1.1 0-1.99.9-1.99 2S5.9 22 7 22s2-.9 2-2-.9-2-2-2zM1 2v2h2l3.6 7.59-1.35 2.45c-.16.28-.25.61-.25.96 0 1.1.9 2 2 2h12v-2H7.42c-.14 0-.25-.11-.25-.25l.03-.12.9-1.63h7.45c.75 0 1.41-.41 1.75-1.03l3.58-6.49c.08-.14.12-.31.12-.48 0-.55-.45-1-1-1H5.21l-.94-2H1zm16 16c-1.1 0-1.99.9-1.99 2s.89 2 1.99 2 2-.9 2-2-.9-2-2-2z"
}), "ShoppingCart");
Fe = le.default = dn;
const un = $()((r) => ({
  searchBox: {
    margin: r.spacing(2)
  },
  menuIcon: {
    marginRight: r.spacing(1),
    marginBottom: 0
  }
})), fn = R(function({ model: e }) {
  const { classes: t } = un(), { selection: i } = e, { pluginManager: a } = z.getEnv(e), s = z.getSession(e), [o, l] = w(), h = a.evaluateExtensionPoint("TrackSelector-multiTrackMenuItems", [], { session: s });
  return f.createElement(
    f.Fragment,
    null,
    i.length ? f.createElement(
      B,
      { className: t.menuIcon, onClick: (p) => l(p.currentTarget) },
      f.createElement(
        ct,
        { badgeContent: i.length, color: "primary" },
        f.createElement(Fe, null)
      )
    ) : null,
    f.createElement(K, { anchorEl: o, open: !!o, onMenuItemClick: (p, _) => {
      _(), l(void 0);
    }, onClose: () => l(void 0), menuItems: [
      { label: "Clear", onClick: () => e.clearSelection() },
      ...h.map((p) => ({
        ...p,
        ..."onClick" in p ? { onClick: () => p.onClick(e) } : {}
      }))
    ] })
  );
}), hn = P(() => import("./FacetedDialog-a4bb0bc1.js")), mn = $()((r) => ({
  searchBox: {
    margin: r.spacing(2)
  },
  menuIcon: {
    marginRight: r.spacing(1),
    marginBottom: 0
  }
}));
function pn({ model: r, setHeaderHeight: e, setAssemblyIdx: t }) {
  const { classes: i } = mn(), [a, s] = w(!1), { filterText: o } = r;
  return f.createElement(
    "div",
    { ref: (l) => e((l == null ? void 0 : l.getBoundingClientRect().height) || 0), "data-testid": "hierarchical_track_selector" },
    f.createElement(
      "div",
      { style: { display: "flex" } },
      f.createElement(sn, { model: r, setAssemblyIdx: t }),
      f.createElement(fn, { model: r }),
      f.createElement(dt, { className: i.searchBox, label: "Filter tracks", value: o, onChange: (l) => r.setFilterText(l.target.value), fullWidth: !0, InputProps: {
        endAdornment: f.createElement(
          ut,
          { position: "end" },
          f.createElement(
            B,
            { onClick: () => r.clearFilterText() },
            f.createElement(ht, null)
          )
        )
      } }),
      f.createElement(ft, { className: i.menuIcon, onClick: () => s(!0) }, "Open faceted selector")
    ),
    f.createElement(xe, { fallback: f.createElement("div", null) }, a ? f.createElement(hn, { handleClose: () => s(!1), model: r }) : null)
  );
}
const gn = R(pn), vn = ({ tree: r, model: e, offset: t }) => typeof jest > "u" ? f.createElement(_t, { disableWidth: !0 }, ({ height: i }) => f.createElement(we, { height: (i || t) - t, model: e, tree: r })) : f.createElement(we, { height: 9e3, model: e, tree: r }), _n = ({ overrideDimensions: r, children: e }) => r ? f.createElement("div", { style: { ...r } }, e) : f.createElement(f.Fragment, null, e), bn = R(function({ model: r, toolbarHeight: e, overrideDimensions: t }) {
  return f.createElement(
    _n,
    { overrideDimensions: t },
    f.createElement(Sn, { model: r, toolbarHeight: e }),
    f.createElement(St, { model: r })
  );
}), Sn = R(function({ model: r, toolbarHeight: e = 0 }) {
  const [t, i] = w(0), [a, s] = w(0), { assemblyNames: o } = r, l = o[t];
  return l ? f.createElement(
    f.Fragment,
    null,
    f.createElement(gn, { model: r, setHeaderHeight: s, setAssemblyIdx: i }),
    f.createElement(vn, { tree: r.hierarchy(l), model: r, offset: e + a })
  ) : null;
}), En = /* @__PURE__ */ Object.freeze(/* @__PURE__ */ Object.defineProperty({
  __proto__: null,
  default: bn
}, Symbol.toStringTag, { value: "Module" }));
export {
  En as H,
  fn as S,
  oe as d
};
