#!/usr/bin/env python3
"""
electrode_viewer.py - versión: toolbar + magnifier por cada figura en Compare

Requisitos: pillow, matplotlib, numpy
"""
import os, io, json
from datetime import datetime
import tkinter as tk
from tkinter import filedialog, messagebox
from PIL import Image, ImageOps
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt

# -------------------------
# Config
# -------------------------
THUMB_SIZE = (220, 140)
DETAIL_THUMB = (220, 120)
DEFAULT_COLUMNS = 5
MAX_SELECT = 3

PHOTO_NAMES = ["thumbnail.png", "thumbnail.jpg", "photo.png", "photo.jpg", "image.png", "image.jpg"]
PULSES_NAMES = ["pulses.png", "pulses.jpg", "pulses.jpeg", "pulses_plot.png", "pulses_plot.jpg"]
AVG_NAMES = ["avg.png", "avg.jpg", "average.png", "average.jpg", "promedio.png", "promedio.jpg"]

# -------------------------
# Helpers
# -------------------------
def find_file_with_any(root, names):
    for n in names:
        p = os.path.join(root, n)
        if os.path.isfile(p):
            return p
    for ext in (".png", ".jpg", ".jpeg", ".bmp"):
        for f in sorted(os.listdir(root)):
            if f.lower().endswith(ext):
                return os.path.join(root, f)
    return None

def parse_date_flexible(s):
    if not s: return None
    if isinstance(s, datetime): return s
    s = s.strip()
    try: return datetime.fromisoformat(s)
    except: pass
    for fm in ("%Y-%m-%d","%d-%m-%Y","%d/%m/%Y","%Y/%m/%d","%d %b %Y","%d %B %Y"):
        try: return datetime.strptime(s, fm)
        except: pass
    import re
    m = re.search(r"(\d{4})[-/\.](\d{1,2})[-/\.](\d{1,2})", s)
    if m:
        y,mo,d = m.groups()
        try: return datetime(int(y), int(mo), int(d))
        except: return None
    return None

def read_metadata(folder):
    out = {"snr": None, "measurement_date": None, "origin_date": None}
    p = os.path.join(folder, "metadata.json")
    if os.path.isfile(p):
        try:
            with open(p, "r", encoding="utf-8") as fh:
                j = json.load(fh)
            if "snr" in j:
                try: out["snr"] = float(j["snr"])
                except: out["snr"] = None
            md = j.get("measurement_date") or j.get("date") or j.get("measured_on")
            od = j.get("origin_date") or j.get("origin") or j.get("created")
            out["measurement_date"] = parse_date_flexible(md) if md else None
            out["origin_date"] = parse_date_flexible(od) if od else None
        except Exception:
            pass
    if out["snr"] is None:
        p2 = os.path.join(folder, "snr.txt")
        if os.path.isfile(p2):
            try:
                with open(p2, "r", encoding="utf-8") as fh:
                    t = fh.read().strip()
                    for token in t.replace(",", " ").split():
                        try:
                            out["snr"] = float(token); break
                        except: continue
            except:
                pass
    return out

def pil_image_to_tk_photoimage(pil_image, size):
    im = pil_image.convert("RGBA")
    im = ImageOps.fit(im, size, Image.LANCZOS)
    b = io.BytesIO(); im.save(b, format="PNG"); b.seek(0)
    return tk.PhotoImage(data=b.read())

# -------------------------
# Model
# -------------------------
class Electrode:
    def __init__(self, path):
        self.path = path
        self.name = os.path.basename(path)
        self.photo_path = find_file_with_any(path, PHOTO_NAMES)
        self.pulses_path = find_file_with_any(path, PULSES_NAMES)
        self.avg_path = find_file_with_any(path, AVG_NAMES)
        meta = read_metadata(path)
        self.snr = meta["snr"]
        self.measurement_date = meta["measurement_date"]
        self.origin_date = meta["origin_date"]

# -------------------------
# Magnifier helper
# -------------------------
class Magnifier:
    def __init__(self, master_tk, fig_canvas, ax, image_path):
        self.master = master_tk
        self.canvas = fig_canvas
        self.ax = ax
        self.image_path = image_path
        self.zoom_factor = 3.0
        self.size_px = 280
        self._connected = False
        self.zoom_win = None
        self._image_arr = None
        self._load_image(image_path)

    def _load_image(self, image_path):
        if image_path and os.path.isfile(image_path):
            try:
                im = Image.open(image_path).convert("RGB")
                self._pil = im
                self._image_arr = np.array(im)
            except:
                self._pil = None
                self._image_arr = None
        else:
            self._pil = None
            self._image_arr = None

    def set_zoom(self, z):
        try:
            v = float(z); self.zoom_factor = max(1.0, v)
        except: pass

    def enable(self):
        if self._image_arr is None:
            messagebox.showinfo("Magnifier", "Imagen no disponible para la lupa.")
            return
        if not self._connected:
            self.cid = self.canvas.mpl_connect("motion_notify_event", self._on_motion)
            self._connected = True
            if self.zoom_win is None or not tk.Toplevel.winfo_exists(self.zoom_win):
                self.zoom_win = tk.Toplevel(self.master)
                self.zoom_win.wm_title("Magnifier")
                self.zoom_label = tk.Label(self.zoom_win)
                self.zoom_label.pack()
                self.zoom_win.geometry(f"{self.size_px}x{self.size_px}")
                self.zoom_win.attributes("-topmost", True)

    def disable(self):
        if self._connected:
            try: self.canvas.mpl_disconnect(self.cid)
            except: pass
            self._connected = False
        if self.zoom_win is not None and tk.Toplevel.winfo_exists(self.zoom_win):
            try: self.zoom_win.destroy()
            except: pass
            self.zoom_win = None

    def update_image(self, image_path):
        self.image_path = image_path
        self._load_image(image_path)

    def _on_motion(self, event):
        if event.inaxes is None or event.inaxes != self.ax: return
        if self._image_arr is None: return
        x = int(event.xdata) if event.xdata is not None else None
        y = int(event.ydata) if event.ydata is not None else None
        if x is None or y is None: return
        h, w = self._image_arr.shape[0], self._image_arr.shape[1]
        half_w = int((self.size_px / self.zoom_factor) / 2)
        half_h = int((self.size_px / self.zoom_factor) / 2)
        x0 = max(0, x - half_w); x1 = min(w, x + half_w)
        y0 = max(0, y - half_h); y1 = min(h, y + half_h)
        crop = self._image_arr[y0:y1, x0:x1]
        if crop.size == 0: return
        pil = Image.fromarray(crop).resize((self.size_px, self.size_px), Image.LANCZOS)
        tkim = pil_image_to_tk_photoimage(pil, (self.size_px, self.size_px))
        if self.zoom_win and tk.Toplevel.winfo_exists(self.zoom_win):
            self.zoom_label.configure(image=tkim); self.zoom_label.image = tkim

# -------------------------
# App
# -------------------------
class App:
    def __init__(self, root):
        self.root = root
        root.title("Electrode Viewer")
        self.base_folder = None
        self.electrodes = []
        self.thumb_cache = {}
        self.selected = {}
        self.columns = DEFAULT_COLUMNS
        self.sort_by = tk.StringVar(value="SNR")

        top = tk.Frame(root); top.pack(fill="x", padx=6, pady=6)
        tk.Button(top, text="Open folder...", command=self.open_folder).pack(side="left")
        tk.Button(top, text="Refresh", command=self.refresh).pack(side="left", padx=4)
        tk.Label(top, text="Columns:").pack(side="left", padx=(12,2))
        self.spin_cols = tk.Spinbox(top, from_=1, to=8, width=3, command=self.on_columns_change)
        self.spin_cols.delete(0,"end"); self.spin_cols.insert(0, str(DEFAULT_COLUMNS)); self.spin_cols.pack(side="left")
        tk.Label(top, text="Sort by:").pack(side="left", padx=(12,2))
        self.sort_menu = tk.OptionMenu(top, self.sort_by, "SNR", "Measurement date"); self.sort_menu.pack(side="left")
        tk.Button(top, text="Apply sort", command=self.build_grid).pack(side="left", padx=6)
        tk.Button(top, text="Compare Selected", command=self.compare_selected).pack(side="right")
        self.info_label = tk.Label(top, text="Seleccione carpeta con subcarpetas de electrodos"); self.info_label.pack(side="left", padx=8)

        self.canvas = tk.Canvas(root, height=700)
        self.v_scroll = tk.Scrollbar(root, orient="vertical", command=self.canvas.yview)
        self.grid_container = tk.Frame(self.canvas)
        self.grid_container.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas.create_window((0,0), window=self.grid_container, anchor="nw")
        self.canvas.configure(yscrollcommand=self.v_scroll.set)
        self.canvas.pack(side="left", fill="both", expand=True)
        self.v_scroll.pack(side="right", fill="y")

    def open_folder(self):
        folder = filedialog.askdirectory(title="Seleccionar carpeta base de electrodos")
        if not folder: return
        self.base_folder = folder; self.load_electrodes(); self.build_grid()

    def refresh(self):
        if not self.base_folder:
            messagebox.showinfo("Info", "Primero abre una carpeta base."); return
        self.load_electrodes(); self.build_grid()

    def load_electrodes(self):
        self.electrodes = []
        if not self.base_folder: return
        for name in sorted(os.listdir(self.base_folder)):
            p = os.path.join(self.base_folder, name)
            if os.path.isdir(p): self.electrodes.append(Electrode(p))
        self.apply_sorting()

    def apply_sorting(self):
        key = self.sort_by.get()
        if key == "SNR":
            self.electrodes.sort(key=lambda x: (-(x.snr if x.snr is not None else -1e9)))
        else:
            def mk(e): return e.measurement_date.timestamp() if e.measurement_date else -1e12
            self.electrodes.sort(key=lambda x: mk(x), reverse=True)

    def on_columns_change(self):
        try: v = int(self.spin_cols.get()); self.columns = max(1, min(8, v))
        except: self.columns = DEFAULT_COLUMNS
        self.build_grid()

    def build_grid(self):
        self.apply_sorting()
        for child in self.grid_container.winfo_children(): child.destroy()
        self.thumb_cache.clear(); self.selected.clear()
        cols = self.columns; pad = 8; r = c = 0
        for e in self.electrodes:
            frame = tk.Frame(self.grid_container, bd=1, relief="groove", padx=6, pady=6)
            frame.grid(row=r, column=c, padx=pad, pady=pad, sticky="n")
            thumb_img = self.make_thumbnail_for_electrode(e, THUMB_SIZE)
            lbl = tk.Label(frame, image=thumb_img, cursor="hand2"); lbl.image = thumb_img; lbl.pack()
            lbl.bind("<Button-1>", lambda ev, ee=e: self.open_detail(ee))
            name_lbl = tk.Label(frame, text=e.name, font=("Arial", 11, "bold")); name_lbl.pack(anchor="w", pady=(6,0))
            snr_text = f"SNR: {e.snr:.2f}" if e.snr is not None else "SNR: N/A"
            md_text = f'Medición: {e.measurement_date.strftime("%Y-%m-%d")}' if e.measurement_date else "Med: N/A"
            od_text = f'Origen: {e.origin_date.strftime("%Y-%m-%d")}' if e.origin_date else "Origin: N/A"
            info_lbl = tk.Label(frame, text=f"{snr_text}  |  {md_text}\n{od_text}", font=("Arial", 9), justify="left"); info_lbl.pack(anchor="w", pady=(4,0))
            var = tk.IntVar(value=0); chk = tk.Checkbutton(frame, text="Select", variable=var); chk.pack(anchor="w", pady=(4,0))
            self.selected[e.name] = (var, e)
            c += 1
            if c >= cols: c = 0; r += 1

    def make_thumbnail_for_electrode(self, e: Electrode, size=THUMB_SIZE):
        if e.path in self.thumb_cache: return self.thumb_cache[e.path]
        if e.photo_path and os.path.isfile(e.photo_path):
            try: pil = Image.open(e.photo_path)
            except: pil = Image.new("RGB", size, (220,220,220))
        else: pil = Image.new("RGB", size, (240,240,240))
        try: tkimg = pil_image_to_tk_photoimage(pil, size)
        except: tkimg = tk.PhotoImage(width=size[0], height=size[1])
        self.thumb_cache[e.path] = tkimg; return tkimg

    # ---------------- Detail window ----------------
    def open_detail(self, electrode: Electrode):
        win = tk.Toplevel(self.root); win.title(f"Detail - {electrode.name}")
        left = tk.Frame(win); left.pack(side="left", fill="y", padx=6, pady=6)
        right = tk.Frame(win); right.pack(side="left", fill="both", expand=True, padx=6, pady=6)
        snr_text = f"SNR: {electrode.snr:.2f}" if electrode.snr is not None else "SNR: N/A"
        md_text = f'Fecha medición: {electrode.measurement_date.strftime("%Y-%m-%d")}' if electrode.measurement_date else "Measurement date: N/A"
        od_text = f'Fecha origen: {electrode.origin_date.strftime("%Y-%m-%d")}' if electrode.origin_date else "Origin date: N/A"
        info_lbl = tk.Label(left, text=f"{electrode.name}\n{snr_text}\n{md_text}\n{od_text}", font=("Arial", 12, "bold"), justify="left"); info_lbl.pack(pady=(8,6))
        thumbs_frame = tk.Frame(left); thumbs_frame.pack(pady=(4,0), fill="y")
        def load_small(path, size):
            if path and os.path.isfile(path):
                try: pil = Image.open(path)
                except: pil = Image.new("RGB", size, (220,220,220))
            else: pil = Image.new("RGB", size, (240,240,240))
            return pil_image_to_tk_photoimage(pil, size)
        photo_thumb = load_small(electrode.photo_path, DETAIL_THUMB)
        pulses_thumb = load_small(electrode.pulses_path, DETAIL_THUMB)
        avg_thumb = load_small(electrode.avg_path, DETAIL_THUMB)
        btn_photo = tk.Button(thumbs_frame, image=photo_thumb, command=lambda: show_on_main(electrode.photo_path, "Photo")); btn_photo.image = photo_thumb; btn_photo.pack(side="top", padx=4, pady=6, anchor="w"); tk.Label(thumbs_frame, text="Photo").pack(side="top", anchor="w")
        btn_pulses = tk.Button(thumbs_frame, image=pulses_thumb, command=lambda: show_on_main(electrode.pulses_path, "Pulses")); btn_pulses.image = pulses_thumb; btn_pulses.pack(side="top", padx=4, pady=6, anchor="w"); tk.Label(thumbs_frame, text="Pulses").pack(side="top", anchor="w")
        btn_avg = tk.Button(thumbs_frame, image=avg_thumb, command=lambda: show_on_main(electrode.avg_path, "Average")); btn_avg.image = avg_thumb; btn_avg.pack(side="top", padx=4, pady=6, anchor="w"); tk.Label(thumbs_frame, text="Average").pack(side="top", anchor="w")

        main_fig = plt.Figure(figsize=(7,6)); main_ax = main_fig.add_subplot(111); main_ax.axis("off")
        main_canvas = FigureCanvasTkAgg(main_fig, master=right); main_canvas.draw(); main_widget = main_canvas.get_tk_widget(); main_widget.pack(fill="both", expand=True)
        toolbar_frame = tk.Frame(right); toolbar_frame.pack(fill="x"); toolbar = NavigationToolbar2Tk(main_canvas, toolbar_frame); toolbar.update()

        # Magnifier controls
        mag_frame = tk.Frame(left); mag_frame.pack(pady=(8,6))
        mag_var = tk.IntVar(value=0)
        mag_check = tk.Checkbutton(mag_frame, text="Magnifier", variable=mag_var); mag_check.pack(side="left", padx=(0,6))
        tk.Label(mag_frame, text="Zoom:").pack(side="left")
        mag_scale = tk.Scale(mag_frame, from_=1.5, to=8.0, resolution=0.5, orient="horizontal", length=140); mag_scale.set(3.0); mag_scale.pack(side="left", padx=(4,0))

        initial_img = electrode.avg_path or electrode.pulses_path or electrode.photo_path
        mag = Magnifier(self.root, main_canvas, main_ax, initial_img)
        def on_mag_toggle(*_):
            if mag_var.get(): mag.set_zoom(mag_scale.get()); mag.enable()
            else: mag.disable()
        mag_check.configure(command=on_mag_toggle); mag_scale.configure(command=lambda v: mag.set_zoom(v))

        def show_on_main(path, title):
            main_ax.clear(); main_ax.axis("off"); main_ax.set_title(title)
            if path and os.path.isfile(path):
                try:
                    im = plt.imread(path); main_ax.imshow(im); mag.update_image(path)
                except Exception:
                    main_ax.text(0.5,0.5,"Cannot load image", ha="center")
            else:
                main_ax.text(0.5,0.5,"No file", ha="center")
            main_canvas.draw()

        if electrode.avg_path: show_on_main(electrode.avg_path, "Average")
        elif electrode.pulses_path: show_on_main(electrode.pulses_path, "Pulses")
        else: show_on_main(electrode.photo_path, "Photo")
        tk.Button(toolbar_frame, text="Reset view", command=lambda: toolbar.home()).pack(side="right", padx=6)

    # ---------------- Compare window ----------------
    def compare_selected(self):
        selected_list = []
        for name,(var,e) in self.selected.items():
            if var.get(): selected_list.append(e)
        if not selected_list:
            messagebox.showinfo("Info", "Seleccioná 1-3 electrodos (casillas) para comparar."); return
        if len(selected_list) > MAX_SELECT:
            messagebox.showwarning("Warning", f"Seleccioná hasta {MAX_SELECT} electrodos."); return
        self.open_compare_window(selected_list)

    def open_compare_window(self, electrodes):
        n = len(electrodes)
        win = tk.Toplevel(self.root); win.title("Compare - " + ", ".join([e.name for e in electrodes]))
        top = tk.Frame(win); top.pack(fill="x", padx=6, pady=6)
        tk.Label(top, text="Show:").pack(side="left")
        view_choice = tk.StringVar(value="All")
        for val in ("Photo","Pulses","Average","All"):
            tk.Radiobutton(top, text=val, variable=view_choice, value=val, indicatoron=0).pack(side="left", padx=4)
        tk.Label(top, text="   Date shown:").pack(side="left", padx=(12,2))
        date_choice = tk.StringVar(value="Medición")
        tk.OptionMenu(top, date_choice, "Medición", "Origen").pack(side="left")
        tk.Button(top, text="Apply", command=lambda: render()).pack(side="right")

        container = tk.Frame(win); container.pack(fill="both", expand=True)

        def date_text(e):
            if date_choice.get() == "Medición":
                return f"Medición: {e.measurement_date.strftime('%Y-%m-%d')}" if e.measurement_date else "Medición: N/A"
            else:
                return f"Origen: {e.origin_date.strftime('%Y-%m-%d')}" if e.origin_date else "Origen: N/A"

        def render():
            for ch in container.winfo_children(): ch.destroy()
            choice = view_choice.get()
            if choice == "All":
                rows = 3; cols = n
                fig = plt.Figure(figsize=(5*cols, 4*rows)); axs=[]
                for r in range(rows):
                    row_axes=[]
                    for c in range(cols):
                        ax = fig.add_subplot(rows, cols, r*cols + c + 1); ax.axis("off"); row_axes.append(ax)
                    axs.append(row_axes)
                for idx,e in enumerate(electrodes):
                    ax_photo = axs[0][idx]
                    if e.photo_path and os.path.isfile(e.photo_path):
                        try: ax_photo.imshow(plt.imread(e.photo_path))
                        except: ax_photo.text(0.5,0.5,"Cannot load photo",ha="center")
                    else: ax_photo.text(0.5,0.5,"No photo",ha="center")
                    snr_txt = f"{e.name}\nSNR: {e.snr:.2f}" if e.snr is not None else f"{e.name}\nSNR: N/A"
                    ax_photo.set_title(f"{snr_txt}\n{date_text(e)}", fontsize=9)
                    ax_p = axs[1][idx]
                    if e.pulses_path and os.path.isfile(e.pulses_path):
                        try: ax_p.imshow(plt.imread(e.pulses_path))
                        except: ax_p.text(0.5,0.5,"Cannot load pulses",ha="center")
                    else: ax_p.text(0.5,0.5,"No pulses file",ha="center")
                    ax_p.set_title("Pulses", fontsize=9)
                    ax_a = axs[2][idx]
                    if e.avg_path and os.path.isfile(e.avg_path):
                        try: ax_a.imshow(plt.imread(e.avg_path))
                        except: ax_a.text(0.5,0.5,"Cannot load avg",ha="center")
                    else: ax_a.text(0.5,0.5,"No avg file",ha="center")
                    ax_a.set_title("Average", fontsize=9)
                canvas = FigureCanvasTkAgg(fig, master=container); canvas.draw(); canvas.get_tk_widget().pack(fill="both", expand=True)
            else:
                # For each electrode, create its own Figure + Canvas + Toolbar + Magnifier controls
                # arrange them horizontally in container
                col_frame = tk.Frame(container); col_frame.pack(fill="both", expand=True)
                for idx,e in enumerate(electrodes):
                    cell = tk.Frame(col_frame, bd=1, relief="flat")
                    cell.pack(side="left", fill="both", expand=True, padx=4, pady=4)

                    fig = plt.Figure(figsize=(6,6))
                    ax = fig.add_subplot(111); ax.axis("off")
                    # load image path depending on choice
                    path = None
                    if choice == "Photo": path = e.photo_path
                    elif choice == "Pulses": path = e.pulses_path
                    elif choice == "Average": path = e.avg_path

                    if path and os.path.isfile(path):
                        try: ax.imshow(plt.imread(path))
                        except: ax.text(0.5,0.5,"Cannot load",ha="center")
                    else:
                        ax.text(0.5,0.5,"No file",ha="center")
                    title = f"{e.name}\nSNR: {e.snr:.2f}" if e.snr is not None else f"{e.name}\nSNR: N/A"
                    title = title + "\n" + date_text(e)
                    ax.set_title(title, fontsize=10)

                    # canvas and toolbar inside this cell
                    canvas = FigureCanvasTkAgg(fig, master=cell)
                    canvas.draw()
                    widget = canvas.get_tk_widget(); widget.pack(fill="both", expand=True)
                    tbf = tk.Frame(cell); tbf.pack(fill="x")
                    toolbar = NavigationToolbar2Tk(canvas, tbf); toolbar.update()

                    # magnifier controls for THIS column (affects this axis)
                    mgf = tk.Frame(cell); mgf.pack(fill="x", pady=(4,4))
                    mag_var = tk.IntVar(value=0)
                    mag_chk = tk.Checkbutton(mgf, text="Magnifier", variable=mag_var)
                    mag_chk.pack(side="left", padx=(4,6))
                    tk.Label(mgf, text="Zoom:").pack(side="left")
                    mag_scale = tk.Scale(mgf, from_=1.5, to=8.0, resolution=0.5, orient="horizontal", length=140)
                    mag_scale.set(3.0); mag_scale.pack(side="left", padx=(4,6))

                    # instantiate magnifier for this axis
                    sample_path = path
                    mag = Magnifier(self.root, canvas, ax, sample_path)
                    def make_toggle(mag_obj, var_obj, scale_obj):
                        def toggle():
                            if var_obj.get():
                                mag_obj.set_zoom(scale_obj.get()); mag_obj.enable()
                            else:
                                mag_obj.disable()
                        return toggle
                    mag_chk.configure(command=make_toggle(mag, mag_var, mag_scale))
                    mag_scale.configure(command=lambda v, m=mag: m.set_zoom(v))

        render()

# -------------------------
def main():
    root = tk.Tk()
    app = App(root)
    root.mainloop()

if __name__ == "__main__":
    main()
