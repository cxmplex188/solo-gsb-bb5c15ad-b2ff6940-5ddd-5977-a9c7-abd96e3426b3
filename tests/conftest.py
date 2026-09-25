"""测试会话级配置：把 SQLite 持久缓存隔离到临时目录。

应用级缓存单例在 ``import app.cache`` 时按环境变量 ``REACTION_CACHE_DB``
建库。conftest 先于任何测试模块被导入，因此在这里把路径指向会话级临时
目录：既不污染仓库内的持久缓存文件，也保证每次运行测试时缓存都是全新
的（``cached`` 标志与全新缓存一致）。
"""

import os
import tempfile

_CACHE_DIR = tempfile.TemporaryDirectory(prefix="reaction-cache-test-")
os.environ["REACTION_CACHE_DB"] = os.path.join(_CACHE_DIR.name, "test_cache.db")
