from cresset import flare
import tkinter as tk
from PIL import ImageGrab
import tempfile
import requests
import pyperclip
import os
from PySide6.QtGui import QIcon

# ===== 配置 =====
SERVER_URL = "http://172.16.1.185:5000/decimer"

# ===== 截图工具（简化版）=====
def take_screenshot():
    """Take a screenshot of a user-selected area."""
    print("👉 拖拽鼠标选择化学结构区域...")
    root = tk.Tk()
    root.attributes("-alpha", 0.3, "-topmost", True, "-fullscreen", True)
    canvas = tk.Canvas(root, bg="black", cursor="cross")
    canvas.pack(fill="both", expand=True)
    coords = {}

    def press(e):
        coords['start'] = (e.x, e.y)
        coords['rect'] = canvas.create_rectangle(e.x, e.y, e.x, e.y, outline='red', width=2)
    def drag(e):
        x1, y1 = coords['start']
        canvas.coords(coords['rect'], x1, y1, e.x, e.y)
    def release(e):
        x1, y1 = coords['start']
        x2, y2 = e.x, e.y
        root.destroy()
        bbox = (min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2))
        global screenshot_image
        screenshot_image = ImageGrab.grab(bbox=bbox)

    canvas.bind("<ButtonPress-1>", press)
    canvas.bind("<B1-Motion>", drag)
    canvas.bind("<ButtonRelease-1>", release)
    root.mainloop()

    return screenshot_image if 'screenshot_image' in globals() else None

# ===== Flare 扩展入口 =====
@flare.extension
class DecimerExtension:
    def load(self):
        """Called by Flare when the extension is loaded."""
        tab = flare.main_window().ribbon["Extensions"]
        tab.tip_key = "X"
        group = tab["Ligand"]

        control = group.add_button("Decimer", self._run_decimer)
        control.icon = QIcon("decimer")   # 自动找 decimer.png
        control.tooltip = (
            "Use DECIMER to recognize chemical structures from screenshots. "
            "Select an area containing a chemical structure, and DECIMER will "
            "convert it to SMILES, copy to clipboard, and add to Flare ligands table."
        )
        control.tip_key = "D"

        print(f"Loaded {self.__class__.__module__}.{self.__class__.__name__}")

    def _run_decimer(self):
        """Ribbon button callback."""
        screenshot_image = take_screenshot()
        if not screenshot_image:
            print("❌ 未截取图像")
            return

        # 保存临时文件
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            screenshot_image.save(tmp.name)
            tmp_path = tmp.name

        try:
            # 发送 HTTP 请求
            with open(tmp_path, 'rb') as f:
                files = {'image': f}
                print("📤 正在发送图像到服务器...")
                resp = requests.post(SERVER_URL, files=files, timeout=30)

            if resp.status_code == 200:
                result = resp.json()
                smiles = result.get("smiles", "").strip()
                if smiles:
                    pyperclip.copy(smiles)
                    print(f"🎉 SMILES 已复制到剪贴板：{smiles}")
                    
                    # 添加到 Flare 项目
                    try:
                        project = flare.main_window().project
                        if project:
                            # 使用 flare.read_string 读取 SMILES 并添加到项目中
                            ligand = flare.read_string(smiles, 'smi')
                            project.ligands.extend(ligand)
                            print(f"✅ SMILES 已添加到 Flare 项目的配体表格中")
                        else:
                            print("⚠️ 未找到 Flare 项目，请确保已打开项目")
                    except Exception as e:
                        print(f"❌ 添加到 Flare 项目失败: {e}")
                        import traceback
                        traceback.print_exc()
                else:
                    print("⚠️ 服务器返回空 SMILES（可能无法识别）")
            else:
                error = resp.json().get("error", "Unknown error")
                print(f"❌ 服务器错误 ({resp.status_code}): {error}")

        except requests.exceptions.RequestException as e:
            print(f"💥 网络错误: {e}")
        except Exception as e:
            print(f"💥 客户端错误: {e}")
        finally:
            os.unlink(tmp_path)
