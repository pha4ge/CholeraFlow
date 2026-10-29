FROM python:3.11-slim
COPY requirements.txt /tmp/choleraflow-requirements.txt
RUN python3 -m pip install --no-cache-dir -r /tmp/choleraflow-requirements.txt \
    && python3 -c "import jinja2; assert jinja2.__version__ == '3.1.6'"
