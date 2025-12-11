#!/usr/bin/env python
from flask import Flask, request, jsonify
from DECIMER import predict_SMILES
import tempfile
import os

app = Flask(__name__)

@app.route('/decimer', methods=['POST'])
def decimer_predict():
    if 'image' not in request.files:
        return jsonify({"error": "No image file provided"}), 400

    file = request.files['image']
    if file.filename == '':
        return jsonify({"error": "Empty filename"}), 400

    try:
        # 保存临时图像
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            file.save(tmp.name)
            tmp_path = tmp.name

        # 调用 DECIMER
        smiles = predict_SMILES(tmp_path)

        # 清理
        os.unlink(tmp_path)

        return jsonify({"smiles": smiles if smiles else ""})

    except Exception as e:
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    print("🚀 启动 DECIMER REST 服务 @ http://0.0.0.0:5000/decimer")
    app.run(host='0.0.0.0', port=5000, threaded=True, debug=False)
    # 仅允许内网访问（可选：绑定 0.0.0.0 对外，但建议限制 IP）
    #app.run(host='0.0.0.0', port=5000, debug=False)
