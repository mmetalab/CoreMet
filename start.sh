#\!/bin/bash
PORT=${PORT:-8080}
WORKERS=${GUNICORN_WORKERS:-2}
TIMEOUT=${GUNICORN_TIMEOUT:-120}

echo "Checking for existing process on port $PORT..."
PIDS=$(lsof -ti :$PORT 2>/dev/null)
if [ -n "$PIDS" ]; then
    echo "Killing existing process(es): $PIDS"
    echo "$PIDS" | xargs kill -9 2>/dev/null
    sleep 0.5
fi

if [ "$CORMET_ENV" = "production" ]; then
    echo "Starting CoreMet (production) on port $PORT with $WORKERS workers..."
    exec gunicorn "app.main:create_app().server" \
        --bind "0.0.0.0:$PORT" \
        --workers "$WORKERS" \
        --timeout "$TIMEOUT" \
        --access-logfile - \
        --error-logfile - \
        --preload
else
    echo "Starting CoreMet (development) on port $PORT..."
    python run.py
fi
